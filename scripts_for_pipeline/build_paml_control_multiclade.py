#!/usr/bin/env python
##build_paml_control.py
##written 11/11/14 by Groves Dixon
ProgramName = 'build_paml_control.py'
LastUpdated = '11/2/17'
#update, added -leafOnly option, which will allow you to label all terminal branches and ancestral nodes only for a polyphyletic forground group

By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Outputs a Paml control file for a given sequence file and a given tree file.

If necessary it collapses taxa in the tree file so that they match those given in 
the sequence file.

Now contains option to annotate the tree file  (adding "#1" to the branch branch leading to a particular clade 
for getting clade-specific dN/dS rates (see PAML_tree_tutorial_molevol2012 by Nicholas R. Meyerson for more on labeling trees for PAML)
'''

AdditionalProgramInfo = '''
Additional Program Information:
Default parameters for the model arguments are intended for running codeml in pairwise mode.

'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio import Phylo
import re
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = False, dest = 'input', help = 'The the input sequence file')
parser.add_argument('-controlName', required = False, default = "input_name + .cnt", dest = 'controlFileName', help = 'The file name to use for the PAML control file you are making')
parser.add_argument('-o', required = False, dest = "out", default = 'input_name + .codeml', help = 'The name you want for the output file when PAML is run using.')
parser.add_argument('-t', required = True, dest = 'tree', help = 'The name of the tree file')
parser.add_argument('-spp', required = True, dest = 'sppList', help = 'A list of ALL species that could be included.')
parser.add_argument('-runMode', required = False, default = '-2', dest = 'runMode', help = 'The run mode you would like to use. The defulat is -2 (pairwise dN/dS calculation).')
parser.add_argument('-model', required = False, default = '0', dest = 'model', help = 'The assignment for the model. Default is 0 (used for pair-wise comparisons).')
parser.add_argument('-NSsites', required = False, default = '0', dest = 'nssites', help = 'Set the values to use for NSsites. (Default = 0)')
parser.add_argument('-fix_omega', required = False, default = '0', dest = 'fixW', help = 'Set the value for fix_omega. (Default = 0)')
parser.add_argument('-omega', required = False, default = '1', dest = 'W', help = 'Set the fixed or initial value for omega. (Default = 0)')
parser.add_argument("-clean_data", required = False, default = '1', dest = 'cleanData', help = "Set to 0 or 1 for whether you want to implement PAML's data cleaning option. This is supposed to clean up your input sequence. The default is to use it, so set to 0 to turn it off.")
parser.add_argument('-suffix', required = False, default = '', dest = 'suffix', help = 'A string to add to the output control file names so they can be differentiated from other control files for that gene. An example would be "AlternativeModel" or just "Alt"')
parser.add_argument('-treeNoteSuffix', required = False, default = '#1', dest = 'treeSuffix', help = 'A string that is added to terminal branches in the tree to designate the lineage for branch models in paml. (See paml manual and Tree_file_notes.txt for how the notations are added. Default value = "#1"')
parser.add_argument('-clade', required = False, default = 'none', dest = 'clade', help = 'A list of species (a subset of that given for -spp) that represent a clade within the tree. It is OK to include species that are missing from the tree, but the set of species must make up a monophyletic group within the tree. The given species will be used to assign the clade to measure dN/dS rates')
parser.add_argument('-minForeground', required = False, default = 2, dest = 'minimum_foreground', help = 'The minimum number of foreground species required to output a control file.')
parser.add_argument('-minBackground', required = False, default = 2, dest = 'minimum_background', help = 'The minimum number of background species required to output a control file.')
args = parser.parse_args()

#Set up some global variables

#set up the file names
File = args.input
TreeFile = args.tree
SppFile = args.sppList
ControlFileName = args.controlFileName
OutFileName = args.out
CladeFile = args.clade
#change the ControlFileName and TreFile names to defaults if non were given
if ControlFileName == "input_name + .cnt":
    ControlFileName = File.rstrip(".codon") + ".cnt"
if OutFileName == 'input_name + .codeml':
    OutFileName = File.rstrip(".codon") + ".codeml"
TreeOutFileName = File.rstrip(".codon") + ".tree"


#set up parameters for running paml
modelType = str(args.suffix)
RunMode = str(args.runMode)
Model = args.model
NSsites = args.nssites
FixOmega = args.fixW
Omega = args.W
CleanData = str(args.cleanData)

#set up arguments for editing the species tree
SuffixInTree = args.treeSuffix
MinForeground = int(args.minimum_foreground)
MinBackground = int(args.minimum_background)


def read_spp_list(SppFile):
    """read in the species list"""
    with open(SppFile, 'r') as infile:
        sppList = []
        for line in infile:
            sppList.append(line.strip("\n"))
    return sppList


def read_clade_file(CladeFile, absenteeList):
    print("Reading in clade file...")
    with open(CladeFile, 'r') as infile:
        cladeDict = {}
        cladeList = []
        cladeSpeciesList = []
        for line in infile:
            line=line.strip("\n").split("\t")
            print(line)
            spp = line[0]
            clade = line[1]
            #skip species that are not present in the codon file
            if spp in absenteeList:
                continue
            try:
                cladeDict[clade].append(spp)
                cladeSpeciesList.append(spp)
            except KeyError:
                cladeList.append(clade)
                cladeDict[clade] = [spp]
                cladeSpeciesList.append(spp)
    print("Found {} clades:".format(len(cladeList)))
    for c in cladeList:
        print("{}".format(c))
        for s in cladeDict[c]:
            print("\t{}".format(s))
    return cladeDict, cladeList, cladeSpeciesList





def read_codon_file(File, sppList):
    """read the codon file and extract the species that have
    sequences so that you can trim the tree file appropritately"""
    print "\nLooking for the following {} species in the sequence file:".format(len(sppList))
    for i in sppList:
        print i
    #set up lists to store the taxa that are present and absent in sequence file
    inSeqFileList = []
    absenteeList = []
    #read through file and find all the taxa 
    with open(File, 'r') as infile:
        for line in infile:
            line = line.strip("\n")
            if line in sppList:
                inSeqFileList.append(line)
    print "\nFound the following {} species in the sequence file:".format(len(inSeqFileList))
    for i in inSeqFileList:
        print i
    #build absentee list that will be removed from the phylogenetic tree to let codeml run for this sequence
    for i in sppList:
        if i not in inSeqFileList:
            absenteeList.append(i)
    print "\nThe following {} species were not found in the sequence file and will be removed from the tree:".format(len(absenteeList))
    for i in absenteeList:
        print i
    return absenteeList

def test_representation(CladeList, absenteeList):
    """Function to make sure we have enough species in the foreground and 
    background lineages"""
    print("\nChecking to make sure there are enough species in the foreground and background lineages...")
    passing = True
    cladePresentList = []
    backgroundPresentList = []
    failMessage = ''
    for i in CladeList:
        if i not in absenteeList:
            cladePresentList.append(i)
    if len(cladePresentList) < MinForeground:
         passing = False
         failMessage += "Less than {} species in foreground clade\n".format(MinForeground)
    for i in sppList:
        if i not in CladeList and i not in absenteeList:
            backgroundPresentList.append(i)
    if len(backgroundPresentList) < MinBackground:
        passing = False
        failMessage += "Less than {} species in background clade\n".format(MinBackground)
    print("\nSpecies in foreground clade with sequences:")
    for sp in cladePresentList:
        print(sp)
    print("\nSpecies in background clade with sequences:")
    for sp in backgroundPresentList:
        print(sp)
    if passing:
        print("\nFound at least {} species in foreground and at least {} species in backgroud. Continuing to output.".format(MinForeground, MinBackground))
    else:
        print(failMessage)
    return passing, failMessage
    

def trim_tree(absenteeList, TreeFile, cladeDict, cladeList, cladeSpeciesList):
    """Collapse away species from the phylogenetic tree that
    are not found in this sequence file. Output the tree file."""
    print "\nReading the Tree..."
    #parse the tree using Phylo
    tree = Phylo.read(TreeFile, 'newick')
    print "Here is the starting tree:"
    Phylo.draw_ascii(tree)
    terminals = tree.get_terminals()
    print "\nFound the following {} taxa in the tree:".format(len(terminals))
    print terminals
    #prune away taxa that are not included for this sequence file
    for taxon in absenteeList:
        print taxon
        tree.prune(taxon)
    print "\nPruned away these species:"
    print absenteeList
    print "\nHere is the tree with the missing taxa pruned away:\n"
    Phylo.draw_ascii(tree)
    #unless you have a clock, PAML requires that your tree is unrooted, ie has a trifurcation at first node. So do that here
    ROOT = tree.get_nonterminals()[0]
    if ROOT.is_bifurcating() == True:
        firstNode = tree.get_nonterminals()[1]
        tree.collapse(firstNode)
    #add notations to the tree to identify the 'foreground' branches
    #these are assigned to a monophyletic group of species assigned with the argument -clade
    #by default add "#1" to the branch leading to the clade. Change -inc from 'no' to make it inclusive,
    #adding #1 to the branch leading to the clade as well as all terminal branches.
    if Model == "2":
        print "\nAssigning the foreground branches in the tree based on the species given in the clade file..."
        print "These are the foreground clade to be labeled:"
        for c in cladeList:
            print("{}".format(c))
        for s in cladeDict[c]:
            print("\t{}".format(s))

        #check if the base of tree and base of foreground clade are same thing:
        print("Assignment done one clade at a time...")
        baseNode = tree.get_nonterminals()[0]         #get the base node
        for clade in cladeList:
            print("\n{}".format(clade))
            cladeSpecies = cladeDict[clade]
            for s in cladeSpecies:
                print("\t{}".format(s))
            baseForClade = tree.common_ancestor(cladeSpecies)#get the base for the foreground clade
            if baseNode.branch_length == baseForClade.branch_length and baseForClade.branch_length == None:
                markedCladeMembers = [] #set up list to record members of the Cladelist that are marked as foreground
                print("\nNote, clade {} was split by trifurcation...".format(clade))
                print("this one will be ignored")
                continue
            if baseForClade.name in cladeSpecies:
                baseForClade.name = baseForClade.name + "#1"
            else:
                baseForClade.name = "#1"


    #if RunMode is not 2 just output the pruned tree as is
    print "\nOutputting the following revised tree for the species content of the sequence file"
    print "it should have a trifurcation at the base unless you are using a clock\n"
    Phylo.draw_ascii(tree)
    # if tree.rooted == False:
    #     print "The revised tree is an unrooted tree (regardless of how the sketch above looks)"
    # if tree.rooted == True:
    #     print "Hmm, the tree is rooted. This may not be right for PAML input. You should check."
    Phylo.write(tree, TreeOutFileName, "newick")
         
#set up the base control file text for Codeml
Base = '''
seqfile = {}
     treefile = {}
      outfile = {}

        noisy = 0   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = {}  * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

    cleandata = 1   * "I added on 07/07/2004" Mikita Suyama

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        model = {}
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = {}   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = {}   * 1: omega or omega_1 fixed, 0: estimate
        omega = {}   * initial or fixed omega, for codons or codon-transltd AAs

    *fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        *alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       *Malpha = 0   * different alphas for genes
        *ncatG = 4   * # of categories in the dG or AdG models of rates

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
        RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
       method = 0   * 0: simultaneous; 1: one branch at a time
'''.format(File, TreeOutFileName, OutFileName, RunMode, Model, NSsites, FixOmega, Omega)




sppList = read_spp_list(SppFile)
#we only need a cladeList when there is a cladeFile
#this is only necessary when setting up a control file for a branch-sites test


absenteeList = read_codon_file(File, sppList)

cladeDict, cladeList, cladeSpeciesList = read_clade_file(CladeFile, absenteeList)

passing, failMessage = test_representation(cladeSpeciesList, absenteeList)
if passing == False:
    ControlFileName = ControlFileName + "_FAIL"
    with open(ControlFileName, 'w') as out:
        out.write(failMessage)
    exit(failMessage)


#only bother trimming the tree if not runnin pairwise
trim_tree(absenteeList, TreeFile, cladeDict, cladeList, cladeSpeciesList)
#output the control file
with open(ControlFileName, 'w') as out:
    out.write(Base)


                
            
            

