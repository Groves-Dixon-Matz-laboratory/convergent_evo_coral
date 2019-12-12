#!/usr/bin/env python
##output_seqs_step1.py
##last updated 12-12-19
    ##revised read_seqs so it can handle relative paths to fasta files
##Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit
from Bio import SeqIO
import datetime
import os
import ntpath





##############################
###### DEFINE FUNCTIONS ######
##############################

def read_seqs(faList):
    '''read in the fasta files'''
    speciesList = []
    print("\nReading in fasta files...")
    seqDict = {}
    for fa in faList:
        faFile=ntpath.basename(fa)
        print("{}...".format(fa))
        spp = faFile.split('_')[0]
        speciesList.append(spp)
        seqDict[spp] = {}
        fasSeqs = SeqIO.parse(open(fa), 'fasta')
        #iterate through the seqs
        for seq in fasSeqs:
            # seqDict[spp][seq.id] = str(seq.seq).upper()
            seq.description=''
            seqDict[spp][seq.id] = seq
    return(speciesList, seqDict)

def read_fastOrtho(orthoFile):
    '''Read through the fastortho .end output and output the sequences for each orthogroup
    '''
    totalGroups = 0
    totalAssigned = 0
    protWritten = 0
    lowRepCount = 0

    #open the fastOrtho groups file
    with open(orthoFile, 'r') as infile:
        #each line contains one orthologous group
        for line in infile:
            totalGroups += 1
            orthoDat = line.strip("\n").split(":")[1]
            oList = orthoDat.split() #the list of the sequences in this orthologous group (with appended speceis names)
            line = line.strip("\n").split()
            groupName = line[0]
            numGenes = line[1].strip("(")
            numTaxa = int(line[2].strip("genes,"))
            if numTaxa < repCut:
                lowRepCount += 1
                continue

            #output the sequences from this orthogroup
            proteinSeqs = []
            for homolog in oList:
                seqid = homolog.split("(")[0]
                fileName = homolog.split("(")[1].split(")")[0]
                taxonName = fileName.split("_")[0]
                
                #gather the protein sequence data
                pseq = protDict[taxonName][seqid]
                proteinSeqs.append(pseq)
            protFileName = "{}/{}.fasta".format(protDir, groupName)
            with open(protFileName, "w") as handle:
                SeqIO.write(proteinSeqs, handle, "fasta")
                protWritten += 1
    print("\n total orthogroups found in {}".format(totalGroups, orthoFile))
    print("{} groups had total taxa below {} and were skipped".format(lowRepCount, repCut))
    print("{} total protein groups written to {}".format(protWritten, protDir))



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
    parser.add_argument('-orthos', required = True, dest = 'ortho_file', help = 'The table of orthologs. Should be tab delimited and should have species names as the column headings that match the species names included in fasta file names')
    parser.add_argument('-prot', required = True, dest = 'prot_fasta_files', nargs="+", help = 'A glob to the protein fastas, for example *PRO.fas')
    # parser.add_argument('-nuc', required = True, dest = 'nuc_fasta_files', nargs="+", help = 'A glob to the protein fastas, for example *PRO.fas')
    # parser.add_argument('-o', required = True, dest = 'output_name',  help = 'Name for the output files')
    parser.add_argument('-ignore', required = False, dest = 'ignore_spp',  help = 'Comma delmited set of species to ignore')
    parser.add_argument('-cut', required = False, default = 5, dest = 'representation_cutoff',  help = 'Minimum number of unique species that must have ortholog to output')


    #----- PARSE ARGUMENTS -----#
    args = parser.parse_args()
    orthoFile = args.ortho_file
    protFaList = args.prot_fasta_files
    # nucFaList = args.nuc_fasta_files
    ignoreSpp = args.ignore_spp
    repCut = int(args.representation_cutoff)
    # outfileName = args.output_name

    #---- SET UP OUTPUT DIR ----#
    now = datetime.datetime.now()
    dirName = 'Orthologs_{}_{}'.format(now.month, now.day)
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    else:
        num=0
        while os.path.exists(dirName):
            num+=1
            dirName = 'Orthologs_{}_{}_{}'.format(now.month, now.day, num)
        os.makedirs(dirName)
    print("\nResults stored in {}".format(dirName))
    
    protDir = "{}/protein_sequences".format(dirName)
    nucDir = "{}/cds_sequences".format(dirName)
    os.makedirs(protDir)
    os.makedirs(nucDir)


    #------ RUN FUNCTIONS ------#

    protSpeciesList, protDict = read_seqs(protFaList)
    # nucSpeciesList, nucDict = read_seqs(nucFaList)
    read_fastOrtho(orthoFile)


    



