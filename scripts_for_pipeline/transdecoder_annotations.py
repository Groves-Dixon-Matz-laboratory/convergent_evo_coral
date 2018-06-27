#!/usr/bin/env python
##transdecoder_annotations.py
##last updated 5-1-18
##Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit
import numpy as np

debug = False

##############################
###### DEFINE FUNCTIONS ######
##############################

def unique(list1):
    x = np.array(list1)
    return(list(np.unique(x)))

def read_annotations():
    print("\nReading in {} annotation files...".format(len(annotFileList)))
    tlineCount = 0
    spDict = {}
    pfDict = {}
    for fileName in annotFileList:
        print("\t{}".format(fileName))
        spp=fileName.split("_")[0]
        with open(fileName, 'r') as infile:
            for line in infile:
                if line[0:5]=="track":
                    continue
                tlineCount += 1
                line = line.strip("\n").split("\t")
                seqName = line[0]
                oline = line[3]
                alist1 = oline.split(",")
                alist2 = alist1[2:]
                if "~~" not in alist1[0]:
                    exit("Error. ~~ expectations not met")
                if "score" not in alist1[1]:
                    exit("Error. score expectation not met")
                sourceList = [x.split("|")[0] for x in alist2]
                tagList = [x.split("|")[1].split(".")[0] for x in alist2]
                nameList = [x.split("|")[2].split(".")[0] for x in alist2]

                #populate lists of the swissprot hits and pfam hits for this sequence
                spList = []
                pfList = []
                for a in alist2:
                    aSplit=a.split("|")
                    if aSplit[0]=='sp':
                        spTag = aSplit[1] + "|" + aSplit[2]
                        if spTag not in spList:
                            spList.append(spTag)
                        else:
                            continue
                    elif aSplit[1][0:2]=="PF":
                        pfTag = aSplit[1] + "|" + aSplit[0]
                        if pfTag not in pfList:
                            pfList.append(pfTag)
                        else:
                            continue
                    else:
                        print(oline)
                        print(a)
                        exit("Error. Not a sp or a pf")
                if len(spList) < 1:
                    spList.append('none')
                if len(pfList) < 1:
                    pfList.append('none')

                #record the swissprot and pfam hits in the dictionary
                spDict[seqName] = spList
                pfDict[seqName] = pfList

                #print out info for debugging
                if debug:
                    print("--------")
                    print(oline)
                    print
                    for x in range(len(tagList)):
                        print("{}\t{}".format(sourceList[x], tagList[x]))
                    print
                    print("SP:")
                    for x in spList:
                        print(x)
                    print("PF:")
                    for x in pfList:
                        print(x)
    return(spDict, pfDict)



def read_orthos():
    print("\nReading in ortholog file {}...".format(infileName))
    groupList = []
    oSpDict = {}
    oPfDict = {}
    with open(infileName, 'r') as infile:
        for line in infile:
            line=line.strip("\n").split()
            groupName = line[0]
            sppSeq = line[1].split(".p")[0]
            spp=sppSeq.split("_")[0]
            if spp in ignoreList:
                continue
            speciesSP = spDict[sppSeq]
            speciesPF = pfDict[sppSeq]
            try:
                oSpDict[groupName] = unique(oSpDict[groupName]+ speciesSP)
                oPfDict[groupName] = unique(oPfDict[groupName] + speciesPF)
            except KeyError:
                groupList.append(groupName)
                oSpDict[groupName] = speciesSP
                oPfDict[groupName] = speciesPF

    print("Found {} orthologous groups".format(len(groupList)))
    return(groupList, oSpDict, oPfDict)



def output():
    with open(outfileName, 'w') as out:
        out.write("orthoGroup\tswissProt\tpfam")
        written = 0
        for o in groupList:
            sp = oSpDict[o]
            pf = oPfDict[o]
            if "none" in sp and len(sp) > 1:
                sp.remove("none")
            if "none" in pf and len(pf) > 1:
                pf.remove("none")
            out.write("\n{}\t{}\t{}".format(o, ";".join(sp), ";".join(pf)))
            written += 1
    print("{} orthogroup annotations written to {}".format(written, outfileName))



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
    parser.add_argument('-i', required = True, dest = 'input_file', help = 'The the input ortholog file')
    parser.add_argument('-abeds', required = True, dest = 'annotation_bed_files', nargs='+', help = 'The bed files output from transdecoder with the annotations')
    parser.add_argument('-ignore', required = False, dest = 'ignore_spp', nargs='+', help = 'Any species in ortholog file that you want to ignore (the first bit in the _ delimited sequence name)')
    parser.add_argument('-o', required = True, dest = 'output_file', help = 'Outputname')
    
    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    infileName = args.input_file
    annotFileList = args.annotation_bed_files
    ignoreList = args.ignore_spp
    outfileName = args.output_file



    #---- RUN FUNCTIONS ----#
    spDict, pfDict = read_annotations()
    groupList, oSpDict, oPfDict = read_orthos()
    output()
    print("\nDone.\n")

