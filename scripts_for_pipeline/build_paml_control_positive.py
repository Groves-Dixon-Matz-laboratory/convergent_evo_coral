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
parser.add_argument('-inc', required = False, default = 'no', dest = 'inclusive', help = 'designate whether you want your "forground" dN/dS value to include terminal branches or just the branch leading to your clade of interest (Default = "no", type "yes" to be inclusive)')
parser.add_argument('-leafOnly', required = False, default = 'no', dest = 'leaf_only', help = 'Change to true if you want to only label terminal branches (leafs) with #1 (not the common ancestor)')
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
Inclusive = args.inclusive
LeafOnly = args.leaf_only
SuffixInTree = args.treeSuffix
MinForeground = int(args.minimum_foreground)
MinBackground = int(args.minimum_background)


def check_arguments(Inclusive, LeafOnly):
    if LeafOnly == "no":
        if Inclusive == "no":
            print("\nYou picked leaf_only = no and inclusive = no")
            print("This means that only the ancestral node for given clade will be labeled with #1")
            print("Your -clade file should be monophyletic")
        else:
            print("\nYou picked leaf_only = no and inclusive = yes")
            print("This means that the ancestral node for your clade will be labeled with #1 as well as all the terminal branches (leaves)")
            print("Your -clade file should be monophyletic")
    else:
        if Inclusive == "no":
            print("\nYou picked leaf_only = yes and inclusive = no")
            print("This option is not allowed. The program will exit")
        else:
            print("\nYou picked leaf_only = yes and inclusive = yes")
            print("This means that all terminal branches from your -clade file will be labeled with #1, as well as the common ancestor to any monophyletic clades within your -clade file")
            print("This will work with a polyphyletic -clade file")




def read_spp_list(SppFile):
    """read in the species list"""
    with open(SppFile, 'r') as infile:
        sppList = []
        for line in infile:
            sppList.append(line.strip("\n"))
    return sppList

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
    

def trim_tree(absenteeList, TreeFile, Inclusive, LeafOnly):
    """Collapse away species from the phylogenetic tree that
    are not found in this sequence file. Output the tree file."""
    print "\nReading the Tree..."
    if LeafOnly != "no" and Inclusive == "no":
        exit("\n\nError. if using 'LeafOnly' labeling you must set Inclusive to 'yes'")
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
        if CladeList != "none":
            if taxon in CladeList:
                CladeList.remove(taxon)
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
        print "These species make up the forground clade:"
        for spp in CladeList:
            print spp
        #identifying the foreground clade works differently depending on whether there are multiple species or just one
        #deal with the case when there are multiple first
        if len(CladeList) > 1:
            print("\nMultiple members of foreground clade remain")
            #add #1 to the node representing the common ancestor to your clade of interest, identifying it as the foreground lineage for the branch sites model
            #this can go one of two ways. When the tree was made unrooted, a trifulrcation was put at the base.
            #This trifurcation can occur within the clade or interest, or not. If it occurs wthin the clade mark the next nodes or branches up in the
            #tree that include the taxa for the clade of interest. Otherwise simply mark the common ancestor of the clase.
            baseNode = tree.get_nonterminals()[0]         #get the base node
            baseForClade = tree.common_ancestor(CladeList)#get the base for the foreground clade
            
            #check if the base of tree and base of foreground clade are same thing:
            if baseNode.branch_length == baseForClade.branch_length and baseForClade.branch_length == None:
                markedCladeMembers = [] #set up list to record members of the Cladelist that are marked as foreground
                print("\nNote, Foreground Clade was split by trifurcation...")
                print("Placing #1 marks higher in tree...")
                #iterate through the nonterminal nodes and check if they include chunks of the foreground clade
                #if they do mark them with #1
                for nonTerminal in tree.get_nonterminals():
                    isForeground = True
                    if nonTerminal.branch_length == None:
                        continue
                    else:
                        for terminalTaxon in nonTerminal.get_terminals():
                            if terminalTaxon.name not in CladeList:
                                isForeground = False
                            else:
                                markedCladeMembers.append(terminalTaxon.name)
                    if isForeground:
                        if LeafOnly == "no":
                            nonTerminal.name = "#1"
                #it's possible that the foreground clade was split into branches going directly from 
                #the base of the tree. In this case they are not all covered by a #1 yet.
                allIncluded = True
                needsMarking = []
                for cMember in CladeList:
                    if cMember not in markedCladeMembers:
                        needsMarking.append(cMember)
                        allIncluded = False
                if not allIncluded:
                    for leaf in tree.get_terminals():
                        if leaf.name in needsMarking:
                            leaf.name = leaf.name + "#1"
            #if the trifurcation did not split foreground clade simply mark the commom ancestor node
            else:
                print("Trifulrcation did not split forground clade")
                if LeafOnly == "no":
                    tree.common_ancestor(CladeList).name = "#1"
            #if you want the foreground lineage to be inclusive for terminal branches, then add the #1s to the terminal taxa in the clade
            if Inclusive != 'no':
                for leaf in tree.get_terminals():
                    if leaf.name in CladeList:
                        leaf.name = leaf.name + "#1"
        #if there is only one member of the clade list left, then it is the sole representative for the lineage, and should be marked #1
        else:
            print("Only a single member of the foreground remains")
            for leaf in tree.get_terminals():
                if leaf.name in CladeList:
                    leaf.name = leaf.name + "#1"

    #new as of 11-22-17
    #wanted to be able to label ancestral nodes of multple clades
    #This only works when you've picked leafsOnly = True
    if LeafOnly != "no" and Inclusive == 'yes':
        print("\n\nAdding hashes to common ancestors of foreground leaves...")
        #iterate through the internal nodes
        for nonTerminal in tree.get_nonterminals():
            isForeground=True
            #for each node, check the terminal branches
            #if any of them do NOT have the #1, then this node is not ancestral to a monophyletic foreground clade
            for terminalTaxon in nonTerminal.get_terminals():
                if terminalTaxon.name[-2:] != "#1":
                    isForeground=False
            if isForeground==True:
                nonTerminal.name = "#1"
    check_arguments(Inclusive, LeafOnly)
    #end of revision



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
if CladeFile != "none":
    CladeList = read_spp_list(CladeFile)
else:
    CladeList = "none"
absenteeList = read_codon_file(File, sppList)

passing, failMessage = test_representation(CladeList, absenteeList)
if passing == False:
    ControlFileName = ControlFileName + "_FAIL"
    with open(ControlFileName, 'w') as out:
        out.write(failMessage)
    exit(failMessage)

#only bother trimming the tree if not runnin pairwise
trim_tree(absenteeList, TreeFile, Inclusive, LeafOnly)
#output the control file
with open(ControlFileName, 'w') as out:
    out.write(Base)


                
            
            

