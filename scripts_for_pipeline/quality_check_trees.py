#!/usr/bin/env python
##quality_check_trees.py
##last upadated 4-24-18
##Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit
from Bio import Phylo
import numpy as np
import os
import datetime


##############################
###### DEFINE FUNCTIONS ######
##############################


def read_in_groups():
    print("\nReading in groups to check monophyly")
    spp2group = {}
    group2spp = {}
    groupList = []
    with open(checkFile, 'r') as infile:
        for line in infile:
            line=line.strip("\n").split()
            spp=line[0]
            group=line[1]
            try:
                group2spp[line[1]].append(line[0])
            except KeyError:
                groupList.append(line[1])
                group2spp[line[1]] = [line[0]]
            
            try:
                spp2group[spp].append(group)
            except KeyError:
                spp2group[spp] = [group]
    print("{} Groups Found:".format(len(groupList)))
    for i in groupList:
        print("{} : {}".format(i, ",".join(group2spp[i])))
    return(groupList, group2spp, spp2group)



def my_is_monophyletic(tree, groupList):
    groupList.sort()
    isMono=False
    singleMiss = False
    targetSize = len(groupList)
    nameList = groupList
    nonTerminals = tree.get_nonterminals()        #all internal nodes
    terminalList = tree.get_terminals()
    allLeafList = [t.name for t in terminalList]
    misplacedList = [] #record the species that don't belong in the clade

    #it's automatically monophyletic if only one species
    if len(groupList)==1:
        isMono=True

    #otherwise iterate through internal nodes and check for monophyly
    else:
        for nt in nonTerminals:
            nodeTerminals = nt.get_terminals()       
            
            #get the species names for the two sets of leaves going in either directio of this node
            nodeLeafList = [t.name for t in nodeTerminals]
            antiNodeLeafList = []
            for leaf in allLeafList:
                if leaf not in nodeLeafList:
                    antiNodeLeafList.append(leaf)
            nodeLeafList.sort()
            antiNodeLeafList.sort()
            if nodeLeafList == groupList:
                isMono=True
            if antiNodeLeafList == groupList:
                isMono=True
    return(isMono)


def check_tree(tree, treeFile, passingList, failGroupList, failedList, failed, verbose):
    passing = True
    with open(treeFile, 'r') as f:
        newickString = f.readline()

    if verbose:
        print("-----------")
        print(treeFile)
        print(Phylo.draw_ascii(tree))
        print("Newick from file:\n{}".format(newickString))

    #prune away anemones, since these make them fail often and arn't trustworthy
    leafList = tree.get_terminals()
    for l in leafList:
        if l.name in ["Aelegantis", "Apallida", "Nvectensis"]:
            tree.prune(l)
            
    treeGroupDict = {}
    groupNameDict = {}
    treeGroupList = []
    leafList = tree.get_terminals()
    if verbose:
        print("Terminal Taxa and their Group Assignments:")
    for leaf in leafList:
        group = spp2group[leaf.name]
        if verbose:
            print("{} : {}".format(leaf.name, group))
        for g in group:
            try:
                treeGroupDict[g].append(leaf)
                groupNameDict[g].append(leaf.name)
            except KeyError:
                treeGroupList.append(g)
                treeGroupDict[g] = [leaf]
                groupNameDict[g] = [leaf.name]
    if verbose:
        print("\nChecking the groups are monophyletic as expected...")
    for g in treeGroupList:
        groupList = groupNameDict[g]
        check = my_is_monophyletic(tree, groupList)
        if verbose:
            print("-")
            print("{} : {}".format(g, ','.join(groupNameDict[g])))
            print(check)
        if check==False:
            passing = False
            failGroupList.append(g)
    if verbose:
        print("Tree Passes: {}".format(passing))
    if passing:
        passingList.append(treeFile.split("RAxML_bestTree.")[1])
    else:
        failed += 1
        failedList.append(treeFile)
    return(passingList, failGroupList, failedList, failed, passing)


def check_trees():
    passingList = []
    failGroupList = []
    faileList = []
    failed = 0
    for treeFile in treeFileList:
        tree = Phylo.read(treeFile, 'newick')
        passingList, failGroupList, failedList, failed, passing = check_tree(tree, treeFile, passingList, failGroupList, faileList, failed, True)
    with open('passingTree.txt', 'w') as out:
        for p in passingList:
            out.write("{}\n".format(p))
    with open('groupFails.txt', 'w') as out:
        for g in failGroupList:
            out.write("{}\n".format(g))

    print("\nOut of {} total trees".format(len(treeFileList)))
    print("{} failed the quality check".format(failed))
    print("{} passed".format(len(passingList)))
    return(failedList, passingList)


def recheck_failures(failedList):
    print("...\n...\n...\nRe-checking the {} failed trees...".format(len(failedList)))
    passingList = []
    failGroupList = []
    faileList = []
    failed = 0
    rescuedList = []
    storedCulpritList = []
    stillFailedList = []
    for treeFile in failedList:
        with open(treeFile, 'r') as f:
            newickString = f.readline()
        tree = Phylo.read(treeFile, 'newick')
        prunePass = False
        nPrunedPass = 0
        resList = []
        culpritList = []
        setLeafList = tree.get_terminals()
        print(setLeafList)
        for i in range(len(setLeafList)):
            # print("pruning {}...".format(setLeafList[i]))
            prunedTree = Phylo.read(treeFile, 'newick')
            leafList = prunedTree.get_terminals()
            # print("Before pruning:")
            # print(Phylo.draw_ascii(prunedTree))
            prunedTree.prune(leafList[i])
            if leafList[i].name != setLeafList[i].name:
                exit("ERROR")
            # print("After pruning:")
            # print(Phylo.draw_ascii(prunedTree))
            passingList, failGroupList, failedList, failed, passing = check_tree(prunedTree, treeFile, passingList, failGroupList, faileList, failed, False)
            if passing:
                prunePass = True
                nPrunedPass += 1
                culpritList.append(leafList[i].name)
                print("RESCUED TREE!")
        if prunePass:
            rescuedList.append(treeFile)
            storedCulpritList.append(culpritList)
        else:
            stillFailedList.append(treeFile)
        print("-------------")
        print("Failed tree {}:".format(treeFile))
        print("Salvaged by pruning? = {}".format(prunePass))
        print("Total prunings that lead to salvaged trees = {}".format(nPrunedPass))
        print("Culprits:")
        for c in culpritList:
            print(c)
        print("Original Newick from file:")
        print(newickString)
    print("\n\n\nTotal resucued = {}".format(len(rescuedList)))
    for i in range(len(rescuedList)):
        print("{} : {}".format(rescuedList[i], ",".join(storedCulpritList[i])))
    return(rescuedList, storedCulpritList, stillFailedList)





def revise_ortho_calls():
    print("\nRevising ortholog calls...")
    toRemoveDict = {}
    totalBadTaxa = 0
    removedOrthoDict = {}
    failedOrthoList = [x.split("RAxML_bestTree.")[1] for x in stillFailedList]
    for i in range(len(rescuedList)):
        treeFileName=rescuedList[i]
        culpritList = storedCulpritList[i]
        ortho = treeFileName.split("RAxML_bestTree.")[1]
        toRemoveDict[ortho] = culpritList
        totalBadTaxa += len(culpritList)

    orthosRevised = []
    totalRemoved = 0
    anemonesRemoved = 0


    outName = orthologFile.split(".txt")[0] + "_revised.txt"
    with open(orthologFile, 'r') as infile:
        with open(outName, 'w') as out:
            for line in infile:
                lineString=line
                line=line.strip("\n").split()
                o = line[0]
                spp = line[1].split("_")[0]
                #first remove all lines for anemones
                if spp in ["Aelegantis", "Apallida", "Nvectensis"]:
                    anemonesRemoved += 1
                    continue
                if o in failedOrthoList:
                    removedOrthoDict[o]=0
                    continue
                try:
                    toRemove = toRemoveDict[o]
                    spp = line[1].split("_")[0]
                    #if it's not an anemone, check if it is a culprit speices and skip the line if it is
                    if spp in toRemove:
                        totalRemoved += 1
                        if o not in orthosRevised:
                            orthosRevised.append(o)
                        continue
                    #write out species that are not culprits or anemones
                    else:
                        out.write(lineString)
                except KeyError:
                    out.write(lineString)
    print("{} total trees failed completely".format(len(stillFailedList)))
    print("{} orthologs from the failed tree group were removed".format(len(removedOrthoDict.keys())))
    print("Of the orthlogs from the rescue group to be revised:")
    print("Total orthologs to revise = {}".format(len(rescuedList)))
    print("Total orthologs revised = {}".format(len(orthosRevised)))
    print("Total taxa to remove = {}".format(totalBadTaxa))
    print("Total removed = {}".format(totalRemoved))
    print("Also removed {} anemone entries from rescued trees".format(anemonesRemoved))
    if len(stillFailedList) == len(removedOrthoDict.keys()) and len(rescuedList) == len(orthosRevised) and totalBadTaxa == totalRemoved:
        print("\tNice, everything makes sense")
    else:
        print("\tHmm, something's wrong. Check match between ortholog file and the trees")





    
##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Check a tree to ensure that a given set of taxa are monophyletic
    '''
    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-trees', required = True, dest = 'tree_files', nargs="+", help = 'Glob to the orthogroup trees')
    parser.add_argument('-groups', required = True, dest = 'check_groups', help = 'Table with the groups to check for monophyly')
    parser.add_argument('-orthos', required = False, default=False, dest = 'ortholog_input', help = 'The ortholog file these trees were built from. eg singleCopyOrthos.txt. Include this if you want to remove discordant taxa from the groups')
    
    #----- PARSE ARGUMENTS -----#
    args = parser.parse_args()
    treeFileList = args.tree_files
    checkFile = args.check_groups
    orthologFile = args.ortholog_input
    # overwrite = args.overwrite


    #------ RUN FUNCTIONS ------#
    groupList, group2spp, spp2group = read_in_groups()
    failedList, passingList = check_trees()
    rescuedList, storedCulpritList, stillFailedList = recheck_failures(failedList)

    
    #print out results summary
    passPct = round(float(len(passingList)) / float(len(treeFileList)), 2)*100
    failPct = round(float(len(failedList)) / float(len(treeFileList)), 2)*100
    rescuePct = round(float(len(rescuedList)) / float(len(failedList)), 2)*100
    totPass = len(rescuedList) + len(passingList)
    totPassPct = round(float(totPass) / float(len(treeFileList)), 2)*100
    print("Total tree checked = {}".format(len(treeFileList)))
    print("Total passing = {} ({}%)".format(len(passingList), passPct))
    print("Total failed = {} ({}%)".format(len(failedList), failPct))
    print("Total failed that could be rescued by pruning a single leaf = {} ({}% of fails)".format(len(rescuedList), rescuePct))
    print("Total passing including rescues = {} ({}%)".format(totPass, totPassPct))

    if orthologFile:
        revise_ortho_calls()



            


