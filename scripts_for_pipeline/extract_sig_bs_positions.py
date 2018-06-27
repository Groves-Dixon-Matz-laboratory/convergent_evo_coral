#!/usr/bin/env python
##vcf_window_wrapper.py
##written 2/7/18 
##Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit



##############################
###### DEFINE FUNCTIONS ######
##############################

def parse_codeml(infileName):
    # print("parsing infile {}...".format(infileName))
    sitesList = []

    with open(infileName, 'r') as infile:
        ready = False
        record = False
        for line in infile:
            if "Bayes Empirical Bayes (BEB)" in line:
                ready=True
            if not ready:
                continue
            if line.strip("\n") == "Positive sites for foreground lineages Prob(w>1):":
                # print("Found start")
                record = True
                continue
            if record:
                if line == "\n":
                    # print("Done.")
                    break
                else:
                    line=line.strip("\n").split()
                    line[2] = line[2].rstrip("*")
                    sitesList.append([infileName]+line)
                    # print(line)
                    continue
            else:
                continue
    return(sitesList)



##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Parses codeml output files, building a table of the sigs that showed evidence of positive selection.
    (Those flagged under "Bayes Empirical Bayes (BEB) analysis" in the output files)
    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, dest = 'input_files', nargs="+", help = 'Glob to alternative model codeml files')
    parser.add_argument('-o', required = True, dest = 'output_file', help = 'Name for output table')
    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    infileList = args.input_files
    outfileName = args.output_file

    recordedFileCount = 0
    totalSitesCount = 0
    totalParsed = 0
    with open(outfileName, 'w') as out:
        out.write("file\tposition\taa\tposterior")
        for infileName in infileList:
            res = parse_codeml(infileName)
            totalParsed += 1
            if totalParsed % 1000 == 0:
                print("\t{} files parsed".format(totalParsed))
            if len(res) > 0:
                recordedFileCount += 1
                totalSitesCount += len(res)
                for siteList in res:
                    out.write("\n" + "\t".join(siteList))
    print("\n{} total files parsed".format(len(infileList)))
    print("{} files had at least one site to record".format(recordedFileCount))
    print("{} total sites recorded".format(totalSitesCount))







    #---- RUN FUNCTIONS ----#

