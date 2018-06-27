#!/usr/bin/env python
##output_seqs_step2.py
##last updated 4-25-18
##Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit
from Bio import SeqIO
import datetime
import os


##############################
###### DEFINE FUNCTIONS ######
##############################


def read_orthos():
    print("\nReading in ortholog file {}...".format(orthoFile))
    oDict = {}
    orthoList = []
    with open(orthoFile, 'r') as infile:
        for line in infile:
            line = line.strip("\n").split()
            try:
                oDict[line[0]].append(line[1])
            except KeyError:
                oDict[line[0]] = [line[1]]
                orthoList.append(line[0])
    return(orthoList, oDict)


def read_seqs(faList):
    """Read in the fasta files in the list into a dictionary with the
    seq records keyed to the seq IDs (which are in the ortholog file)."""
    print("\nReading in fasta files...")
    seqDict = {}
    sppDict = {}
    for fa in faList:
        print("{}...".format(fa))
        spp = fa.split('_')[0]
        fasSeqs = SeqIO.parse(open(fa), 'fasta')
        #iterate through the seqs
        for seq in fasSeqs:
            # seqDict[spp][seq.id] = str(seq.seq).upper()
            seq.description=''
            seqID = seq.id
            speciesName = seqID.split("_")[0]
            sppDict[speciesName] = 0
            seq.id = speciesName #change the sequence ID to just the species name
            seqDict[seqID] = seq #record in dictionary based on ortholog tag
    sppList = sppDict.keys()
    sppList.sort()
    return(seqDict, sppList)


def write_seqs(seqDict, outDir, sppList):
    print("\nWriting out sequences to {}...".format(outDir))
    written = 0
    failedRep = 0
    tot = 0
    for o in orthoList:
        tot += 1
        outName = "{}/{}.fasta".format(outDir, o)
        oSeqNameList = oDict[o]
        if len(oSeqNameList) >= repCut:
            seqList = []
            for s in oSeqNameList:
                if "XP_" in s:
                    continue
                seqList.append(seqDict[s])
            SeqIO.write(seqList, outName, 'fasta')
            written += 1

        else:
            failedRep += 1
    print("{} total orthologs considered".format(tot))
    print("{} failed representation cutoff of >={} terminal taxa".format(failedRep, repCut))
    print("{} were written".format(written))
    if sum([written, failedRep]) == tot:
        print("All orthologs accounted for.")
    else:
        print("Error. Lost some orthologs!")



##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    '''
    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-orthos', required = True, dest = 'ortho_file', help = '')
    parser.add_argument('-prot', required = True, dest = 'prot_fasta_files', nargs="+", help = 'A glob to the protein fastas, for example *PRO.fas')
    parser.add_argument('-nuc', required = True, dest = 'nuc_fasta_files', nargs="+", help = 'A glob to the protein fastas, for example *PRO.fas')
    parser.add_argument('-rcut', required = True, dest = 'representation_cutoff', help = 'Cutoff for how many species must be present to use the ortholog')
    parser.add_argument('-odir', required = False, default=False, dest = 'output_directory', help = 'Directory to put output in. Default = paml_month_day')


    #----- PARSE ARGUMENTS -----#
    args = parser.parse_args()
    orthoFile = args.ortho_file
    protFaList = args.prot_fasta_files
    nucFaList = args.nuc_fasta_files
    repCut = int(args.representation_cutoff)
    oDir = args.output_directory
    # outfileName = args.output_name

    #---- SET UP OUTPUT DIR ----#
    if oDir:
        dirName = oDir
    else:
        now = datetime.datetime.now()
        dirName = 'paml_{}_{}'.format(now.month, now.day)
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    print("\nResults stored in {}".format(dirName))
    protDir = "{}/protein_sequences".format(dirName)
    nucDir = "{}/cds_sequences".format(dirName)
    codonDir = "{}/codon_sequences".format(dirName)
    if not os.path.exists(protDir):
        os.makedirs(protDir)
    if not os.path.exists(nucDir):
        os.makedirs(nucDir)


    #------ RUN FUNCTIONS ------#
    orthoList, oDict = read_orthos()
    protDict, sppList = read_seqs(protFaList)
    nucDict, sppList = read_seqs(nucFaList)
    write_seqs(protDict, protDir, sppList)
    write_seqs(nucDict, nucDir, sppList)


    



