#!/usr/bin/env python
##concatenate_genes_into_nexus.py
##written 6/26/14 by Groves Dixon
ProgramName = 'concatenate_genes_into_nexus.py'
LastUpdated = '10/23/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script takes a set of .codon files output and concatenates them into
a nexus file (.nex). 

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
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
# parser.add_argument('-i', required = False, dest = 'input', help = 'The the input file')
parser.add_argument('-output', '-o', required = True, dest = 'out', help = 'The desired name for the output file')
parser.add_argument('-f', required = True, dest = 'files', nargs = "+", help = 'A glob to the codon files output from pal2nal')
parser.add_argument('-spp', required = True, dest = 'sppList', help = 'A text file listing all the species included.')
parser.add_argument('-l', required = False, default = 10, dest = 'lengthCut', help = 'The length an alignment must be in order to be included')
args = parser.parse_args()

#Assign Arguments
FileList = args.files
OutfileName = args.out
SpeciesFile = args.sppList
LengthCut = int(args.lengthCut)

def read_sppList(SpeciesFile):
    """Read in the species List"""
    speciesList = []
    with open(SpeciesFile, 'r') as infile:
        for line in infile:
            speciesList.append(line.strip('\n'))
    print "\nBuilding Concatenated Sequences for the following Species:"
    for i in speciesList:
        print i
    return speciesList
            
def filter_fileList(FileList):
    for file in FileList:
        with open(file, 'r') as infile:
            lineNumber = 0
            speciesCount = 0
            for line in infile:
                line = line.strip('\n')
                lineNumber += 1
                if lineNumber == 1:
                    line = line.split()
                    #get the count for species in this file
                    sppCount = int(line[0])
                    if sppCount == 0:
                        print "\nSkipping File {} Because it Has no Data".format(file)
                        FileList.remove(file)

def read_files(FileList, speciesList):
    """Read though each of the files and build up a 
    sequence dictionary with each species's sequence 
    linked to the the species name."""
    seqDict = {}
    charSetDict = {}
    totalLength = 0 #to record the total length of the concatenated genes
    lengthFail = 0
    writeOutList = []
    for spp in speciesList:
        seqDict[spp] = '' #set up an empty string for each species
    for file in FileList:
        #set up a list of the species that have not been accounted for in this file (initially all of them)
        #these will be removed incrementally as sequenes for the species included in the file are read in
        #then missing data will be plugged in for all species in the species list not found in the file
        sppLeft = []
        for i in speciesList:
            sppLeft.append(i)
        passLength = True
        # print sppLeft
        # print
        # print file
        with open(file, 'r') as infile:
            lineNumber = 0
            speciesCount = 0
            for line in infile:
                line = line.strip('\n')
                lineNumber += 1
                if lineNumber == 1:
                    line = line.split()
                    #get the count for species in this file
                    sppCount = int(line[0])
                    #get the sequence length
                    seqLength = int(line[1])
                    if seqLength < LengthCut:
                        lengthFail += 1
                        passLength = False
                        continue
                    #record the location of this gene in the concatenated sequence in a dictionary
                    charSetDict[file] = "{}-{}".format((totalLength + 1), (totalLength + seqLength))
                    #set the new totalLength
                    totalLength = totalLength + seqLength
                    # print sppCount
                    # print seqLength
                    continue
                elif line.split("_")[0] in speciesList:
                    #when we hit a new species, start adding sequence for it
                    speciesCount += 1
                    activeSpp = line.split("_")[0]
                    sppLeft.remove(activeSpp)
                    continue
                else:
                    #append the new sequence onto the concatenated sequence for this species
                    seqDict[activeSpp] += line
        if passLength:
            #append missing data for species that did not have a gene sequence from this particular file
            for unincluded in sppLeft:
                seqDict[unincluded] += 'N'*seqLength
            ##check that all sequences are the same length
            for i in speciesList:
                if len(seqDict[i]) != totalLength:
                    print "Warning, species {} doesn't have the right length".format(i)
            writeOutList.append(file)
    print("{} total sequences considered".format(len(FileList)))
    print("{} of these failed the length cutoff of {}".format(lengthFail, LengthCut))
    return seqDict, charSetDict, totalLength, writeOutList

def output(OutfileName, HeadString, speciesList, seqDict, charSetDict, totalLength, writeOutList):
    """Output the results as a nexus file"""
    with open(OutfileName, 'w') as out:  
        out.write(HeadString)
        for spp in speciesList:
            datString = "\n{}\t{}".format(spp, seqDict[spp])
            out.write(datString)
        out.write("\n;\nend;\nbegin assumptions;")
        for i in writeOutList:
            out.write("\n\tcharset {} = {};".format(i, charSetDict[i]))   
        out.write("\nend;")
         
                    
            
        


speciesList = read_sppList(SpeciesFile)
filter_fileList(FileList)
seqDict, charSetDict, totalLength, writeOutList = read_files(FileList, speciesList)

HeadString = """\
#NEXUS
begin data;
dimensions ntax={} nchar={};
format datatype=dna interleave=no missing=N gap=-;
matrix""".format(len(speciesList), totalLength)

output(OutfileName, HeadString, speciesList, seqDict, charSetDict, totalLength, writeOutList)
    

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


