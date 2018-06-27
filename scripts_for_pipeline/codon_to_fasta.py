#!/usr/bin/env python

#Get_files
'''
'''
import math
from sys import argv

debug = False

#make list of files that fit argument
fileName = argv[1]

with open(fileName, 'r') as infile:
    lineNumber = 0
    sequenceLineCount = 0
    for line in infile:
        lineNumber += 1
        if lineNumber == 1:
            seqLength = int(line.strip("\n").split()[1])
            seqLineThreshold = int(math.ceil(float(seqLength) / 60.0))
            continue
        else:
            line = line.strip("\n")
            #this means the line is a title line
            if sequenceLineCount == 0:
                print(">" + line)
                sequenceLineCount += 1
                continue
            elif sequenceLineCount < seqLineThreshold:
                print(line)
                sequenceLineCount += 1
                continue
            elif sequenceLineCount == seqLineThreshold:
                print(line)
                sequenceLineCount = 0
                continue
            else:
                exit('Error!')
            
                
        

