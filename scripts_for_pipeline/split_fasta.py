#!/usr/bin/env python
##split_fasta.py
##written 6/26/14 by Groves Dixon
ProgramName = 'split_fasta.py'
LastUpdated = '6/5/16'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Takes a fasta input file and splits it into a given 
number of equally sized sub fastas, then writes a 
command file for reciprocally blasting them.

'''


##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
from Bio import SeqIO
import numpy as np
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = False, dest = 'input', help = 'The the input file')
parser.add_argument('-n', required = False, dest = 'splits', help = 'The the number of times to split the file')
args = parser.parse_args()

#Assign Arguments
infileName = args.input
splits = int(args.splits)

def split_fasta(infileName):
    '''Function to read in a file as a list of lists
    '''
    count = 0
    fasSeqs = SeqIO.parse(open(infileName), 'fasta')
    for seq in fasSeqs:
        count += 1
    seqsPer = count / splits
    seqSets = range(1, (count - seqsPer), seqsPer)
    seqSets.append((count+1))
    print("\nFound {} fasta sequences in file".format(count))
    print("\nSplitting these seqs over files like so:")
    print(seqSets)
    seqCount = 0
    setIndex = 1
    prefix = fileName = ''.join(infileName.split(".")[0:-1])
    outFileName = prefix + "_split" + str(setIndex) + ".fasta"
    out = open(outFileName, 'w')
    switchFilesAt = seqSets[setIndex]
    fasSeqs = SeqIO.parse(open(infileName), 'fasta')
    for seq in fasSeqs:
        seqCount += 1
        if seqCount == 1:
            out.write(">" + seq.id + "\n" + str(seq.seq))
            continue
        if seqCount < switchFilesAt:
            out.write("\n>" + seq.id + "\n" + str(seq.seq))
        else:
            print("Finished Writing file {}".format(outFileName))
            out.close()
            setIndex += 1
            try:
                switchFilesAt = seqSets[setIndex]
                outFileName = outFileName = prefix + "_split" + str(setIndex) + ".fasta"
                out = open(outFileName, 'w')
                out.write("\n>" + seq.id + "\n" + str(seq.seq))
            except IndexError:
                break
                


split_fasta(infileName)

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


