#!/usr/bin/env python
##paraPruner.py
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

def read_lengths():
    ldict = {}
    with open(lengthsInput, 'r') as infile:
        for line in infile:
            line=line.strip("\n").split()
            ldict[line[0]] = float(line[1])
    return ldict

def setup_dirs():
    now = datetime.datetime.now()
    mainDir = "pruned_results_{}_{}".format(now.month, now.day)
    num = 0
    if overwrite == False:
        while os.path.exists(mainDir):
    		num += 1
    		mainDir = "pruned_results_{}_{}_{}".format(now.month, now.day, num)
    singleCopyDir = "{}/single_copy".format(mainDir)
    cSingleCopyDir = "{}/collapsed_single_copy".format(mainDir)
    subtreeDir =  "{}/subtree".format(mainDir)
    if not os.path.exists(mainDir):
        os.makedirs(mainDir)
    if not os.path.exists(singleCopyDir):
        os.makedirs(singleCopyDir)
    print("\nResults will be stored in {}".format(mainDir))
    return(mainDir, singleCopyDir)

def unique(list1):
    x = np.array(list1)
    return(list(np.unique(x)))

def read_in_trees():
    print("\nReading in trees...")
    emptyFileList = []
    treeList = []
    for treeFile in treeFileList:
        try:
            tree = Phylo.read(treeFile, 'newick')
        except ValueError:
            emptyFileList.append(treeFile)
            print("Warning, tree file {} doesn't seem to have a tree.".format(treeFile))
            continue
        treeList.append(tree)
        tree.name = treeFile.split(".newick")[0]
    print("\t{} tree files were empty, suggest checking these:".format(len(emptyFileList)))
    for t in emptyFileList:
        print(t)
    print("\t{} total trees for analysis".format(len(treeList)))
    return treeList


def find_single_copy(treeList):
    singleCopyList = []
    hasParaList = []
    for tree in treeList:
        terminalSeqList = [term.name for term in tree.get_terminals()]
        terminalSppList = [term.name.split("_")[0] for term in tree.get_terminals()]
        uSppList = unique(terminalSppList)
        singleCopy = len(uSppList) == len(terminalSppList)
        if singleCopy:
            singleCopyList.append(tree)
        else:
            hasParaList.append(tree)
        # print("-----------")
        # print(tree.name)
        # print(terminalSeqList)
        # print(terminalSppList)
        # print(uSppList)
        # print("Single Copy = {}".format(singleCopy))
    print("\t{} total trees checked for paralogs".format(len(treeList)))
    print("\t{} trees were single copy".format(len(singleCopyList)))
    print("\t{} got paralogs".format(len(hasParaList)))
    return(singleCopyList, hasParaList)


def get_longest_seq(leafList):
    """use sequence lengths from read_lengths() to pick longest
    version for an single-species clade"""
    nameList = [leaf.name for leaf in leafList]     #get the leaf names, which are the sequence names
    lengthList = [ldict[name] for name in nameList] #get the sequence lengths pulled with read_lengths()
    for i in reversed(range(len(lengthList))):      #reversed so first becomes assigned if they are equal
        if lengthList[i] == max(lengthList):
            longestSeq = leafList[i]
    # print("------------")
    # print(leafList)
    # print(nameList)
    # print(lengthList)
    # print("longest = {}".format(longestSeq.name))
    return(longestSeq)

def collapse_monophyletic(treeList):
    print("\nLooking for single-species clades...")
    collapseCount = 0
    collapsedList = []
    for tree in treeList:
        nodeList = tree.get_nonterminals()
        for node in nodeList:
            nodeTerminals = node.get_terminals()
            nodeSpp = [leaf.name.split("_")[0] for leaf in nodeTerminals]
            uNodeTaxa = unique(nodeSpp)
            # print('--')
            # print("N terminals = {}".format(len(nodeTerminals)))
            # print("N unique = {}".format(len(uNodeTaxa)))

            #use the unique set of taxa in this node as a test for a single-species clade
            if len(uNodeTaxa) == 1:
                # print("===============================")
                # print("Found one!")
                # Phylo.draw_ascii(tree)
                nameList = [leaf.name for leaf in nodeTerminals]
                toKeep = get_longest_seq(nodeTerminals)
                toRemove = nodeTerminals
                toRemove.remove(toKeep)
                for leaf in toRemove:
                    tree.prune(leaf)
                collapseCount += 1
                # print('------------------------')
                # Phylo.draw_ascii(tree)

        collapsedList.append(tree)
    print("\t{} single-species clades were collapsed into single sequence".format(collapseCount))
    return(collapsedList)



def prune_uninportant(treeList):
    print("\nPruning away anemone species and Adig references...")
    unimportant = ["Aelegantis", "Apallida", "Nvectensis", "XP"]
    prunedUnimport = 0
    prunedPolytomy = 0
    revisedTrees = []
    for tree in treeList:
        originalTree=tree
        leafList = tree.get_terminals()
        for leaf in leafList:
            if leaf.name.split("_")[0] in unimportant:
                tree.prune(leaf)
                prunedUnimport += 1
        revisedTrees.append(tree)

    #removal of these branches can result in polytomies
    #these can then have multiple instances of single species in them
    #so remove those as if they were single-species clades as in collapse_monophyletic()
    for tree in treeList:
        nodeList = tree.get_nonterminals()
        for node in nodeList:
            if node.is_preterminal():
                leafList = node.get_terminals()
                leafSppList = [leaf.name.split("_")[0] for leaf in leafList]
                uSppList = unique(leafSppList)

                #check if any species appears more than once in the polytomy
                if len(leafSppList) > len(uSppList):
                    # print("FOUND POLYTOMOY WITH REPEATED SPECIES")
                    # Phylo.draw_ascii(tree)

                    sppCountList = [leafSppList.count(spp) for spp in leafSppList]
                    repeatSppLeafs = []
                    for i in range(len(sppCountList)):
                        if sppCountList[i] > 1:
                            repeatSppLeafs.append(leafList[i])

                    toKeep = get_longest_seq(repeatSppLeafs)
                    toRemove = repeatSppLeafs
                    toRemove.remove(toKeep)
                    # print("LEAFS TO REMOVE:")
                    # print(uSppList)
                    # print(leafSppList)
                    # print(sppCountList)
                    # print(repeatSppLeafs)
                    # print(toRemove)
                    for leaf in toRemove:
                        tree.prune(leaf)
                        prunedPolytomy += 1
                    # print("---------------")
                    # Phylo.draw_ascii(tree)
                else:
                    continue
            else:
                continue
    print("\t{} leafs from unimportant group pruned away".format(prunedUnimport))
    print("\t{} repeated species within resulting polytomies were pruned".format(prunedPolytomy))
    return(revisedTrees)



def remove_poorly_supported(treeList):
    print("\nPruning away repeated species branches that are not well supported...")
    for tree in treeList:
        nodeList = tree.get_nonterminals()
        print(nodeList)
        for node in nodeList:
            if node.is_preterminal():
                leafList = node.get_terminals()
                leafSppList = [leaf.name.split("_")[0] for leaf in leafList]
                uSppList = unique(leafSppList)
                print(node)





def output_revised_groups(treeList):
    outName = 'singleCopyOrthos.txt'
    print("\nSingle copy orthos saved as {}...".format(outName))
    with open(outName, 'w') as out:
        for tree in treeList:
            leafList = tree.get_terminals()
            seqNameList = [leaf.name for leaf in leafList]
            for sn in seqNameList:
                out.write("{}\t{}\n".format(tree.name, sn))
        



##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Look through gene trees made from alignments of sequences in orthologous groups.
    Output trees and sequence sets for single-copy orthologs.
    Then collapse clades that are monophyletic for a given species and keep only the longest sequence to get more.
    '''
    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-trees', required = True, dest = 'tree_files', nargs="+", help = 'Glob to the orthogroup trees')
    parser.add_argument('-l', required = True, dest = 'length_file', help = 'A table of lengths for all the OFRs used in the blast')
    
    #----- PARSE ARGUMENTS -----#
    args = parser.parse_args()
    treeFileList = args.tree_files
    lengthsInput = args.length_file
    # overwrite = args.overwrite


    #------ RUN FUNCTIONS ------#
    treeList0 = read_in_trees()

    # remove_single_paralogs(treeList0)
    # exit()


    #read in trees and get intial single-copy orthogroups
    print("\nGathering single copy trees...")
    rawSingleCopy, remainingTreeList = find_single_copy(treeList0)
    
    #collapse single-species clades, keeping only longest sequence
    #(asseme that these are isoforms that failed to cluster in cd-hit)
    ldict = read_lengths()
    collapsedTreeList = collapse_monophyletic(remainingTreeList)
    print("\nGathering additional single copy trees after collapsing single-species nodes...")
    collapseSingleCopy, remainingTreeList = find_single_copy(collapsedTreeList)
    allSingleCopy = rawSingleCopy + collapseSingleCopy
    print("\t{} total single copy orthogroups".format(len(allSingleCopy)))
    

    #prune unimportant branches
    prunedTrees = prune_uninportant(remainingTreeList)
    print("\nGathering additional single copy trees after pruning anemones and Adig reference...")
    postPruneSingleCopy, remainingTreeList = find_single_copy(prunedTrees)
    allSingleCopy = rawSingleCopy + collapseSingleCopy + postPruneSingleCopy
    print("\t{} total single copy orthogroups".format(len(allSingleCopy)))


    #prune repeated taxa if they are nearly in-paralogs
    #(ie its bifurcation from itself is poorly supported)
    # remove_poorly_supported(remainingTreeList)

    

    output_revised_groups(allSingleCopy)


    print("\nWriting out remainging trees...")
    if not os.path.exists('with_paralogs'):
        os.makedirs('with_paralogs')
    for r in remainingTreeList:
        Phylo.write(r, "with_paralogs/{}_postPrune.newick".format(r.name), 'newick')


            


