#!/usr/bin/env python
##filter_blast_for_orthos.py 
##Groves Dixon
##last revised 4-24-18

#import modules
import time
import argparse
from sys import argv
from sys import exit




##############################
###### DEFINE FUNCTIONS ######
##############################

def read_lengths():
    ldir = {}
    with open(lengthsInput, 'r') as infile:
        for line in infile:
            line=line.strip("\n").split()
            ldir[line[0]] = float(line[1])
    return ldir

def filter():
    lineCount = 0
    filtered = 0
    written = 0
    with open(outfileName, 'w') as out:
        with open(infileName, 'r') as infile:
            for line in infile:
                lineCount += 1
                lineList = line.split()
                alignLength = float(lineList[3])
                sLength = ldir[lineList[1]]
                if (alignLength/sLength) >= cutoff:
                    out.write(line)
                    written += 1
                else:
                    filtered += 1
                if lineCount %1000000 == 0:
                    print("{} lines processed".format(lineCount))
    print("\n{} total lines processed".format(lineCount))
    print("{} hit lenghts were less than {}% of the subject lengths and were removed".format(filtered, round(cutoff, 3)*100))
    print("{} hits were retained".format(written))
    if filtered + written == lineCount:
        print("All hits accounted for.")
    else:
        print("Error. Lost some hits!")
    print("\nDone.")




##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Use the input ORF lengths to filter a set of blast results
    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, dest = 'input_file', help = 'The the blast results (assumes default -outfmt 6')
    parser.add_argument('-l', required = True, dest = 'length_file', help = 'A table of lengths for all the OFRs used in the blast')
    parser.add_argument('-c', required = False, default=0.75, dest = 'coverage_cut', help = 'Float giving the coverage cutoff to keep a blast hit (alignmentLength / subjectLength must be greater than this). default=0.75')
    parser.add_argument('-o', required = True, dest = 'output_file', help = 'Name for output')

    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    infileName = args.input_file
    lengthsInput = args.length_file
    cutoff = float(args.coverage_cut)
    outfileName = args.output_file

    #---- RUN FUNCTIONS ----#
    ldir = read_lengths()
    filter()







