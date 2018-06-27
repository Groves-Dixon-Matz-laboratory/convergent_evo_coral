#!/usr/bin/env python
## parse_rst.py
##written 8/1/16 by Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio import Phylo
import os
debug = False

#Functions

def read_clade_file(cladeFile):
    """Function to read in the taxa you want to test for convergent evolution between.
    This shoudl be a two-column table with the species in column 1 and the clade name
    in column 2.
    The output dictionaries will be used for all the included .rst files.
    """
    print("\n\nReading clade file {}...".format(cladeFile))
    speciesList = []
    cladeNameList = []
    spp2cladeDict = {}
    clade2sppDict = {}
    with open(cladeFile, 'r') as infile:
        for line in infile:
            line=line.strip("\n").split()
            spp = line[0]
            cladeName = line[1]
            speciesList.append(spp)
            spp2cladeDict[spp] = cladeName
            try:
                clade2sppDict[cladeName].append(spp)
            except KeyError:
                clade2sppDict[cladeName] = [spp]
                cladeNameList.append(cladeName)
    print("\nFound {} species in file:".format(len(speciesList)))
    print("divided among {} clades that have the trait of interest:".format(len(cladeNameList)))
    for c in cladeNameList:
        print("\t{}:".format(c))
        print("\t\t{}".format(" ".join(clade2sppDict[c])))

    return speciesList, cladeNameList, spp2cladeDict, clade2sppDict

def gather_changes(infileName):
    '''Function to parse the rst file and gather the amino acid change data.
    These are recorded as lists keyed to the number for the daughter node of the branch.



    EG for branch number 3 in Rod Page's Tree,
    with a Y>F substitution at position 125 that occured in that branch with 98% posterior:
        {3: ['125', 'Y', '0.980', 'F']}

    '''
    print("parsing file {}...".format(infileName))

    #varialbes for location in file
    foundTree = 0
    recordedTree = 0
    changesSection = 0
    currentBranchNumber = 0
    branchesDone = 0

    #dictionaries to store amino acid changes along branches
    #these are based on the Rod Page's TreeView tree, with each branch numbered
    changesDict = {}       #stores the AA changes found along each branch as a list keyed to branch number
    sppDict = {}           #stores the species each branch is associated with, or 'internal' if its not a terminal branch
    nodesDict = {}         #stores the nodes that are connected by each branch
    daughterNodeList = []


    #read through file
    with open(infileName, 'r') as infile:

        #first thing is tree. Record that and otherwise skip lines
        for line in infile:
            if "Nodes" in line and "are ancestral" in line:
                rootNode = line.split()[1]

            if "tree with node labels" in line:
                foundTree = 1
                continue
            elif foundTree == 1:
                rodsTree = line
                print("\nFound tree with node lables:\n" + rodsTree)
                foundTree = 0
                recordedTree=1
                treeName = infileName + "rodTree.newick"
                with open(treeName, 'w') as out:
                    out.write(line)
                rodTreeString = line
                continue
            elif "Summary of changes along branches" in line:
                if recordedTree == 0:
                    exit("Error. Did not find tree with node labels")
                else:
                    changesSection = 1

            #after tree are the branch summary data
            #these end at a line List of extant and reconstructed sequences, so stop recording there

            #the following records amino acid change data 

            elif "Branch" in line and branchesDone == 0:
                previousBranch = currentBranchNumber
                currentBranchNumber = int(line.split(":")[0].split("Branch ")[1])
                nodes = line.split()[2].split("..")
                daughterNode = nodes[1]
                if daughterNode not in daughterNodeList:
                    daughterNodeList.append(daughterNode)
                else:
                    exit("ERROR. BRANCH LEADING TO MORE THAN ONE DAUGHTER NODE")



                if currentBranchNumber != 1+previousBranch:
                    exit("Error. Branches have not come in proper order")

                print("\nFound new branch number = {}:".format(currentBranchNumber))
                print("Between nodes {} and {}".format(nodes[0], nodes[1]))
                print("Will work with this branch based on the number of its daugher node: {}".format(daughterNode))
                currentBranch = int(daughterNode)
                if len(line.split()) > 3:
                    spp = line.split()[3].split("(")[1].split(")")[0]
                    print("Leading to taxon {}".format(spp))
                else:
                    spp='internal'
                changesDict[currentBranch] = []
                sppDict[currentBranch] = spp
                nodesDict[currentBranch] = nodes
            elif currentBranchNumber != 0 and branchesDone == 0:

                #This is signal line that branch data are complete
                if "List of extant and reconstructed sequences" in line:
                    print("\nFinished Collecting Summary of Changes on Branches data.")
                    currentBranch = 0
                    branchesDone = 1
                    continue
                else:
                    line=line.strip("\n").split()
                    if len(line)<1:
                        continue
                    else:
                        pos=line[0]
                        oldAA = line[1]
                        prob = line[2]
                        newAA = line[4]
                        changeDat = [pos, oldAA, prob, newAA]
                        changesDict[currentBranch].append(changeDat)
                    continue
            elif recordedTree == 0:
                continue
            else:
                continue

    branchList = changesDict.keys()
    branchList.sort()
    print("\nDone parsing file.")
    print("Found {} branches in tree".format(len(branchList)))
    for b in branchList:
        print("\nBranch {}:".format(b))
        print("Nodes = {}".format(nodesDict[b]))
        print("Spp = {}".format(sppDict[b]))
        for d in changesDict[b]:
            print(d)


    return treeName, rootNode, changesDict, nodesDict, sppDict, rodTreeString




def tree(treeName, rootNode, targetSpeciesList, spp2cladeDict, cladeNameList, rodTreeString, clade2sppDict):
    """Read in the Rod Page's TreeView tree and assign the branches to the clades
    given in the clade file. This will allow assignment of the Amino acid changes
    to the particular clades with the phenotype of interest, which are given in the 
    .rst file attached to branch numbers, as given in the Rod Pages TreeView tree.

    Outputs are clade2branchDict and branch2cladeDict, which have the branches linked to the 
    clades of interest in each direction.


    """
    #Read in the Rod Page's tree
    print("\nAsessing tree...")
    tree = Phylo.read(treeName, 'newick')
    originalTree = tree
    print "Here is the starting tree:"
    Phylo.draw_ascii(tree)
    terminalList = tree.get_terminals()
    nonTerminalList = tree.get_nonterminals()
    sppInTreeList = [str(leaf).split("_")[1] for leaf in terminalList]

    #check the root of the tree against expectation from rst file
    root=nonTerminalList[0].confidence
    print("root node in tree here is {}".format(root))
    print("root based on first ancestral node given in rst file is {}".format(rootNode))
    if str(root) == str(rootNode):
        print("Good they match. Directionality of amino acid changes will assume this root")
        print("(note that you can't see the node labels in the printout here)")
    else:
        exit("ERROR. The root found in tree does not match expectation")


    #Assign the termainal branches (leaves) to clades based on the species names
    print("\nAssigning each terminal branch to a clades...")
    clade2leafDict = {}       #record the list of leaf names from tree linked to each target clade
    branch2cladeDict = {}     #record which clade a branch is assigned to. Use branch number as given in .rst file
    clade2branchDict = {}     #record all the terminal and internal branches linked to the clades
    treeTargetList = []       #record target species from clade file that are present in tree
    treeCladeList = []        #record target clades from clade file that are present in tree
    num2nameDict = {}         #record the terminal and internal branch names keyed to the species or clade name (species for leaves, clade for internal)
    for leaf in terminalList:
        leafNum = int(str(leaf).split("_")[0])
        leafSpp = str(leaf).split("_")[1]
        num2nameDict[leafNum] = leafSpp
        if leafSpp in targetSpeciesList:
            treeTargetList.append(leafSpp)
            leafClade = spp2cladeDict[leafSpp] #assign this leaf's clade based off it's species name
            leafCladeSppList = clade2sppDict[leafClade] #grab the list of species for this particular clade

            #next we want to label this as the 'All' clade branch if 
            #the leaf is the only representative for its clade
            inTreeCladeSppList = []
            for cladeSpp in leafCladeSppList:
                if cladeSpp in sppInTreeList:
                    inTreeCladeSppList.append(cladeSpp)
            if len(inTreeCladeSppList) == 1:
                num2nameDict[leafNum] = "All_{}_single_".format(leafClade) + leafSpp #rewrite the branch name indicating this is the only species for this clade
        else:
            leafClade = "background"
        branch2cladeDict[leafNum] = leafClade
        try:
            clade2leafDict[leafClade].append(leaf.name)
            clade2branchDict[leafClade].append(leafNum)
        except KeyError:
            clade2leafDict[leafClade] = [leaf.name]
            clade2branchDict[leafClade] = [leafNum]
            treeCladeList.append(leafClade)
        print("\t{}: {}".format(leaf, leafClade))


    #gather target clades found in tree
    treeTargetCladeList = []
    for t in treeCladeList:
        if t != "background":
            treeTargetCladeList.append(t)

    #Assign the internal nodes to clades if they are monophyletic for them
    print("\nAssigning internal branches to clades...")
    for node in nonTerminalList:
        targetNode = False
        nodeAssignment = "background"
        internalBranchNumber = int(node.confidence)
        nodeLeafList = node.get_terminals()           
        nodeNameList = [n.name for n in nodeLeafList] #record leaf names
        nodeNameList.sort()                           #sort for comparison
        
        if debug:
            print("--------")
            print(node)
            print(internalBranchNumber)
            print(nodeLeafList)
            print(leafNameList)

        #check if clade contains only members of one of the clades of interest
        #ie all species in the node name list fall within one of the clade name lists
        for clade in treeTargetCladeList:
            cladeNameList = clade2leafDict[clade] #recall names for the clade
            cladeNameList.sort()

            cladeSpecific = True
            for sppName in nodeNameList:
                if sppName not in cladeNameList:
                    cladeSpecific = False
            if cladeSpecific:
                print("\tnode {} contains only species from clade {}".format(internalBranchNumber, clade))
                targetNode = True
                if nodeNameList == cladeNameList:
                    print("it's monophletic for the clade")
                    containedSppString = "All_" + clade
                else:
                    print("it's paraphyletic for the clade")
                    nodeSppList = [sp.split("_")[1] for sp in nodeNameList]
                    containedSppString = "-".join(nodeSppList)
                nodeAssignment = containedSppString + "_lineage"
                branch2cladeDict[internalBranchNumber] = clade
                clade2branchDict[clade].append(internalBranchNumber)
        if not targetNode:
            if debug: print("node is background")
        num2nameDict[internalBranchNumber] = nodeAssignment


    #report summary
    print("\nDone assigning branches to clades")
    print("{} species from tree come from clades of interest:".format(len(treeTargetList)))
    for s in treeTargetList: print("\t"+s)
    print("divided accross {} clades of interest present in tree.".format(len(treeTargetCladeList)))
    print(treeName)
    print("Rod Page's TreeView:")
    print(rodTreeString)
    print("Final Branch Number assignments:")
    for c in treeTargetCladeList:
        print("\t{} branches:".format(c))
        for b in clade2branchDict[c]:
            bName = num2nameDict[b]
            print("\t\t{} - {}".format(b, bName))


    return treeTargetCladeList, treeTargetList, clade2branchDict, branch2cladeDict, num2nameDict



def make_sub_dict(changeList):
    positionList = [x[0] for x in changeList]
    changeDict = {}
    for j in range(len(positionList)):
        pos = positionList[j]
        subDat = changeList[j]
        changeDict[pos] = subDat
    return(positionList, changeDict)



def find_convergence(treeTargetCladeList, changesDict, clade2branchDict, branch2cladeDict, num2nameDict, geneId, out):
    """Function to search through the changes that occured along tree and identify overlapping substitutions"""
    posPairs = possible_pairs(treeTargetCladeList)
    print("\nSeaching for convergence events...")
    print("\t{} clades were found in tree".format(len(treeTargetCladeList)))
    print("\t{} possible pairing(s)".format(posPairs))
    pairsCompared = 0
    for i in range(len(treeTargetCladeList)):
        otherCladeList = treeTargetCladeList[(i+1):]
        clade1 = treeTargetCladeList[i]
        clade1Branches = clade2branchDict[clade1]
        for clade2 in otherCladeList:
            pairsCompared += 1
            clade2Branches = clade2branchDict[clade2]
            print("\tComparing clade {} and clade {}...".format(clade1, clade2))
            print("\tComparison {} of {} for this tree.".format(pairsCompared, posPairs))
            for b1 in clade1Branches:
                for b2 in clade2Branches:
                    b1Name = num2nameDict[b1]
                    b2Name = num2nameDict[b2]
                    b1ChangeList = changesDict[b1]
                    b2ChangeList = changesDict[b2]
                    b1PositionList, b1ChangeDict = make_sub_dict(b1ChangeList)
                    b2PositionList, b2ChangeDict = make_sub_dict(b2ChangeList)
                    for pos in b1PositionList:
                        if pos in b2PositionList:
                            b1SubDat = [str(b1), b1Name, branch2cladeDict[b1]] + b1ChangeDict[pos]
                            b2SubDat = [str(b2), b2Name, branch2cladeDict[b2]] + b2ChangeDict[pos]
                            resList = [geneId]
                            for z in range(len(b1SubDat)):
                                resList.append(b1SubDat[z])
                                resList.append(b2SubDat[z])
                            print(resList)
                            out.write("\n" + "\t".join(resList))
    print("Done searching for convergence.")


def factorial(n):
    f=1
    for i in range(1, n+1):
        f=f*i
    return(f)

def possible_pairs(cladeKeys):
    N = len(cladeKeys)
    r = 2
    res = factorial(N) / (factorial(N-r) * factorial(r))
    return(res)



def ref_clade_pair_text(rawCladePair):
    cladePair0 = []
    for r in rawCladePair:
        newR = []
        for r2 in r:
            newR.append(r2.split("_")[1])
        cladePair0.append(".".join(newR))
    cladePair = "__".join(cladePair0)
    return cladePair



#main
def main():
    Start_time = time.time() ##keeps track of how long the script takes to run
    

    ##Set Up Argument Parsing
    Description = '''
    Description: This script parses a .rst file output from codeml ancestral state reconstruction.
    It is assumed paml was run on a control file made from build_paml_ancestralRecon.py
    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, dest = 'input_list', nargs="+", help = 'The the input rst file (generated by codeml)')
    parser.add_argument('-c', required = True, dest = "clade_file", help = "A text file of clades to look for convergent evolution between. One line per clade. Multiple terminal branches can be included as comma separated list. Each clade should be monophyletic. The branch leading to each clade's common ancestor is assumed to be an opportunity for convergent evolution relative to the other clades.")    # parser.add_argument('-o', required = True, dest = 'out', help = 'The desired name for the output file')
    parser.add_argument('-o', required = True, dest = "output_name", help = "Desired name for output file")
    args = parser.parse_args()


    inputList = args.input_list
    cladeFile = args.clade_file
    outfileName = args.output_name

    with open(outfileName, 'w') as out:
        resultsHeader = "ortholog\tb1Num\tb2Num\tb1Name\tb2Name\tb1Clade\tb2Clade\tb1pos\tb2pos\tb1Anc\tb2Anc\tb1post\tb2post\tb1New\tb2New"
        out.write(resultsHeader)
        targetSpeciesList, cladeNameList, spp2cladeDict, clade2sppDict = read_clade_file(cladeFile)
        print("\nParsing data for {} input .rst files...".format(len(inputList)))
        for i in range(len(inputList)):
            infileName = inputList[i]
            geneId = infileName.split(".rst")[0]
            print("\n\n----------------------------------")
            print("Reading input file {} ... (file {} of {})".format(infileName, i+1, len(inputList)))
            treeName, rootNode, changesDict, nodesDict, sppDict, rodTreeString = gather_changes(infileName)           
            treeTargetCladeList, treeTargetList, clade2branchDict, branch2cladeDict, num2nameDict = tree(treeName, rootNode, targetSpeciesList, spp2cladeDict, cladeNameList, rodTreeString, clade2sppDict)
            find_convergence(treeTargetCladeList, changesDict, clade2branchDict, branch2cladeDict, num2nameDict, geneId, out)

    
    #return time to run
    Time = time.time() - Start_time
    print('\nTime took to run: {}'.format(Time))        


if __name__ == "__main__":
   main()   


