#!/usr/bin/env python
##parse_codeml_dnds.py
##written 2/3/15 by Groves Dixon
ProgramName = 'parse_codeml_dnds.py'
LastUpdated = '2/3/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Pulls the second number f
w (dN/dS) for branches:  0.03034 0.13653


'''

AdditionalProgramInfo = '''
Additional Program Information:


'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = False, dest = 'input', help = 'The the input file')
parser.add_argument('-output', '-o', required = True, dest = 'out', help = 'The desired name for the output file')
args = parser.parse_args()

#Assign Arguments
InfileName = args.input
OutfileName = args.out

def read_file():
    '''Function to read in a file as a list of lists
    '''
    lineNumber = 0
    geneList = []
    npList = []
    likelihoodList = []
    with open(InfileName, 'r') as infile:
        for line in infile:
            lineNumber += 1
            line = line.strip('\n').split()
            try:
                npList.append(line[4].strip('):'))
                likelihoodList.append(line[5])
                geneList.append(line[0])
            except IndexError:
                print("Failed to parse the following line:")
                print(line)
                continue
    with open(OutfileName, 'w') as out:
        out.write("contig\tnumParameters\tlikelihood")
        for i in range(len(geneList)):
            out.write("\n{}\t{}\t{}".format(geneList[i], npList[i], likelihoodList[i]))

read_file()

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


