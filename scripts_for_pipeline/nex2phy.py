#!/usr/bin/env python
##nex2phy.py
##written 6/26/14 by Groves Dixon

Description = '''
Description:
Uses biopython to convert an alignment in nexus format to phylip format
'''

AdditionalProgramInfo = '''
Additional Program Information:
'''

##Import Modules 

import argparse
from sys import argv
from sys import exit
from Bio import AlignIO

##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = True, dest = 'input', help = 'The nexus file')
parser.add_argument('-o', required = False, default = "none_given", dest = 'out', help = 'The desired name for the output file name. Default = inputName.phy')
parser.add_argument('-partitions', required = False, default = "y", dest = 'part', help = 'Default is to output a partitions file along with the phylip file. Set this argument to "n" to disable this feature.')
args = parser.parse_args()

#Assign Arguments
InfileName = args.input
OutfileName = args.out
Part = args.part
if OutfileName == "none_given":
    OutfileName = "".join(InfileName.split(".")[0:-1]) + ".phy"



def convert(InfileName, OutfileName):
    '''Uses Biophython to convert the file
    '''
    count = AlignIO.convert(InfileName, "nexus", OutfileName, "phylip")
    print("\nConverted %i alignments" % count)
    print("\nOutput saved as %s" % OutfileName)

def handle_partitions(InfileName, OutfileName):
    """Gets the partitions from the input nexus file and 
    outputs a partitions file for to go along with the phylip output"""
    partitionOutName = "".join(OutfileName.split(".")[0:-1]) + "_partitions.txt"
    partitionList = []
    with open(InfileName, 'r') as infile:
        record = False
        for line in infile:
            if line == 'begin assumptions;\n':
                record = True
                continue
            if record == False:
                continue
            if 'charset' in line:
                partition = line.split()[-1]
                partition = partition.strip(";\n")
                partitionList.append(partition)
    with open(partitionOutName, 'w') as out:
        count = 0
        for i in partitionList:
            count += 1
            if count == 1:
                out.write("DNA, p{}={}".format(count, i))
            else:
                out.write("\nDNA, p{}={}".format(count, i))
    print("\nFound %i partitions in the nexus file" % count)
    print("\nOutput partitions file as %s" %partitionOutName)
    
            
def main():
    print "\n\nConverting nexus file {} to phylip file {}...".format(InfileName, OutfileName)
    convert(InfileName, OutfileName)
    if Part == "y":
        handle_partitions(InfileName, OutfileName)
    print
main()


