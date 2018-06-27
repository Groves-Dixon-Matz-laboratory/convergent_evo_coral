#!/usr/bin/env python
##print_tree.py
##written 3/13/15 by Groves Dixon
ProgramName = 'print_tree.py'
LastUpdated = '3/13/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Reads in a tree using biopython and prints it to standard output.
Intended to run on a lot of trees to scan through them manually.
'''

AdditionalProgramInfo = '''
Additional Program Information:

'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio import Phylo
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = True, dest = 'input', nargs = "+", help = 'The the input tree file(s)')
args = parser.parse_args()


Infiles= args.input

def trim_tree(Infiles):
    """reads in a tree file and prints it"""
    #parse the tree using Phylo
    for treefile in Infiles:
        print "---------------------------"
        print "Tree File = {}".format(treefile)
        tree = Phylo.read(treefile, 'newick')
        print tree
        Phylo.draw_ascii(tree)

trim_tree(Infiles)

                
            
            

