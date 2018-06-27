#!/usr/bin/env python
##build_paml_control.py
##written 11/11/14 by Groves Dixon
ProgramName = 'build_paml_control.py'
LastUpdated = '1/27/15'
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
parser.add_argument('-i', required = True, dest = 'input', help = 'The the input sequence file')
parser.add_argument('-controlName', required = False, default = "input_name + .cnt", dest = 'controlFileName', help = 'The file name to use for the PAML control file you are making')
parser.add_argument('-o', required = True, dest = "out", default = 'input_name + .codeml', help = 'The name you want for the output file when PAML is run using.')
parser.add_argument('-t', required = True, dest = 'tree', help = 'The name of the tree file')
parser.add_argument('-spp', required = True, dest = 'sppList', default="speciesList.txt", help = 'A list of ALL species that could be included.')
args = parser.parse_args()


#assign arguments
File = args.input
TreeFile = args.tree
SppFile = args.sppList
ControlFileName = args.controlFileName
OutFileName = args.out
#change the ControlFileName to the default if none was given
if ControlFileName == "input_name + .cnt":
    ControlFileName = File.split(".rst")[0] + ".cnt" #######for some reason I do not understand using strip here removes the c from beginning of Amil contig names
#change the ControlFileName to the default if none was given
if OutFileName == 'input_name + .codeml':
    OutFileName = File.split(".rst")[0] + ".mcl"
TreeOutFileName = File.split(".rst")[0] + ".tree"


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
            lineSpp = line.split("_")[0]
            if lineSpp in sppList:
                inSeqFileList.append(lineSpp)
    print "\nFound the following {} species in the sequence file:".format(len(inSeqFileList))
    for i in inSeqFileList:
        print i
    if len(inSeqFileList) < 2:
        print("\n\nError. The sequence file {} has only {} species in it.".format(File, len(inSeqFileList)))
        exit("Since there are not enough terminal taxa to maintain the phylogeny exiting")
    #build absentee list that will be removed from the phylogenetic tree to let codeml run for this sequence
    for i in sppList:
        if i not in inSeqFileList:
            absenteeList.append(i)
    print "\nThe following {} species were not found in the sequence file and will be removed from the tree:".format(len(absenteeList))
    for i in absenteeList:
        print i
    return absenteeList

def trim_tree(absenteeList, TreeFile):
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
        print("Removing absent")
        tree.prune(taxon)
    print "\nPruned away these species:"
    print absenteeList
    print "\nHere is the tree with the missing taxa pruned away:\n"
    Phylo.draw_ascii(tree)


    #unless you have a clock, PAML requires that your tree is unrooted, ie has a trifurcation at first node. So do that here
    # ROOT = tree.get_nonterminals()[0]
    # if ROOT.is_bifurcating() == True:
    #     firstNode = tree.get_nonterminals()[1]
    #     tree.collapse(firstNode)
    
    
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
seqfile = {}    * sequence data filename
treefile = {}   * tree structure file name
outfile = {}    * main result file name

noisy = 9    * 0,1,2,3,9: how much rubbish on the screen
verbose = 2  * 0: concise; 1: detailed, 2: too much
runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
             * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 3  * 1:codons; 2:AAs; 3:codons-->AAs

clock =  0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
aaRatefile = /home1/02260/grovesd/bin/paml4.8/dat/jones.dat  * only used for aa seqs with model=empirical(_F)

       * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

model = 2
           * models for codons:
           * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
           * models for AAs or codon-translated AAs:
           * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
           * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

icode = 4  * 0:universal code; 1:mammalian mt; 2-10:see below
Mgene = 0     * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff             
       * AA: 0:rates, 1:separate

fix_alpha = 0   * 0: estimate gamma shape parameter; 1: fix it at alpha
alpha = 0.5     * initial or fixed alpha, 0:infinity (constant rate)
Malpha = 1      * different alphas for genes
ncatG = 4       * # of categories in dG of NSsites models


getSE = 0         * 0: don't want them, 1: want S.E.s of estimates

RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

Small_Diff = .5e-6
cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
method = 1     * Optimization method 0: simultaneous; 1: one branch a time
'''.format(File, TreeOutFileName, OutFileName)

sppList = read_spp_list(SppFile)
#we only need a cladeList when there is a cladeFile
#this is only necessary when setting up a control file for a branch-sites test



absenteeList = read_codon_file(File, sppList)
#only bother trimming the tree if not runnin pairwise
trim_tree(absenteeList, TreeFile)
#output the control file
with open(ControlFileName, 'w') as out:
    out.write(Base)


                
            
            

