#!/usr/bin/env python
##longest_isotig.py
##written 11/13/17 by Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit
# import numpy as np
from Bio import SeqIO


#Functions
def read_cdhit(infileName):
    print("\n\n\nReading in cd-hit data {}".format(infileName))
    '''Function to read in the cd-hit data into a dictionary
    '''
    clusterDict = {}
    nCluster = 0
    nLongest = 0
    nShort = 0
    with open(infileName, 'r') as infile:
        for line in infile:
            lineList = line.strip('\n').split()
            if lineList[0] == ">Cluster":
                nCluster += 1
                currentCluster = lineList[1]
                continue
            elif "*" in line:
                contig = line.split("... *")[0].split(">")[1]
                try:
                    clusterDict[contig]
                    print("-----------------------------")
                    print(line)
                    print(clusterDict[contig])
                    print(currentCluster)
                    print(contig)
                    exit("Error! Repeat contig names. Note the cd-hit can't handle seq ids longer than 19 characters")
                except KeyError:
                    clusterDict[contig] = currentCluster
                    nLongest += 1
            else:
                nShort += 1
                continue
    if nCluster != len(clusterDict.keys()):
        print("N clusters read = {}".format(nCluster))
        print("Number of keys assigned in dictionary = {}".format(len(clusterDict.keys())))
        exit("\nHmmm. The dictionary key length is not of expected length. Check inputs.")
    else:
        print("\nFound {} clusters in cd-hit input file".format(nCluster))
        print("Found {} longest sequences".format(nLongest))
        print("Recorded {} longest sequences in the dictionary (should be equal)".format(len(clusterDict.keys())))
        print("Will not write out {} shorter isotigs".format(nShort))
        print("{} + {} = {} total sequences".format(nShort, nLongest, sum([nShort, nLongest])))
    return clusterDict, nShort


def output_seqs(infileName, outfileName, clusterDict, nShort):
    '''Use the clusterDictionary to read through fasta
    writing out only the longest isotig from each cluster'''
    written = 0
    skipped = 0
    print("Writing out results to file {}...".format(outfileName))
    fasSeqs = SeqIO.parse(open(infileName), 'fasta')
    #iterate through the seqs
    with open(outfileName, 'w') as out:
        for seq in fasSeqs:
            subName =str(seq.id)[0:19]
            try:
                clusterNumber = clusterDict[subName]
                out.write(">{}\n".format(seq.id))
                out.write(str(seq.seq) + "\n")
                written += 1
            except KeyError:
                skipped += 1
                # print("SKIP")
                # print(seq.id)
                continue
    print("\nWrote {} longest isotig sequences".format(written))
    print("Skipped {} short isotigs".format(skipped))
    if written == len(clusterDict.keys()):
        print("Good, wrote out one sequence for each cluster")
        print(" Note empty sequence idenities are also skipped, so may have additional skipped sequences.")
    else:
        print("ERROR!, number of written sequences does not match expectation. Check inputs")
    # if skipped == nShort:
    #     print("Good, skipped expected number of sequences")
    # else:
    #     print("ERROR!, number of skipped sequences does not match expectation. Check inputs")





#main
def main():
    Start_time = time.time() ##keeps track of how long the script takes to run
    

    ##Set Up Argument Parsing
    Description = '''
    Description:
    Read through output from cd-hit and writes only the top contig for each isogroup to a new fasta file
    These are indcated by an asterisk in the cd-hit output.
    Eg:
    cd-hit-est -i transcriptome.fasta -o transcriptome_clust.fasta -c 0.99 -G 0 -aL 0.3 -aS 0.3
    longest_isotig.py -i transcriptome.fasta -cdh transcriptome_clust.fasta -o transcriptome_longest.fa
    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = False, dest = 'infile_name', help = 'The the input transcriptome')
    parser.add_argument('-cdh', required = False, dest = 'cdhit_results', help = 'The cd-hit results')
    parser.add_argument('-o', required = True, dest = 'outfile_name', help = 'The desired name for the output file')
    args = parser.parse_args()


    infileName = args.infile_name
    cdhitInput = args.cdhit_results
    outfileName = args.outfile_name

    

    clusterDict, nShort = read_cdhit(cdhitInput)
    output_seqs(infileName, outfileName, clusterDict, nShort)



    
    #return time to run
    Time = time.time() - Start_time
    print('\nTime took to run: {}'.format(Time))        



if __name__ == "__main__":
   main()   


