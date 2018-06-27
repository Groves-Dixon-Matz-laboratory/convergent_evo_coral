#!/usr/bin/env python

#gather_go_terms.py
#loops through the idmappings for uniprot
#collects all GO terms for all uniprot annotations 
#associated with a GO term in a table

from sys import argv

debug = False

#set up usage information
usage = "gather_go_terms.py singleCopyAnnotations.tsv idmapping_selected.tab\n"

try:
    argv[1]
except IndexError:
    print('\nUsage:')
    exit(usage)




infileName = argv[1]
mappings = argv[2]


#read in the ortholog annotations
#make a dictionary with the uniprot annotations keyed to lists of 
#all associated orthologs
print("\nReading in ortholog annotations...")
oList = []
aList = []
oDict = {}
aDict = {}
notAnnotated=0
totalOrtho = 0
with open(infileName, 'r') as infile:
	for line in infile:
		totalOrtho += 1
		line=line.strip("\n").split("\t")
		o = line[0]
		spList = line[1].split(";")
		annotList = [a.split("|")[0] for a in spList]
		aList += annotList
		if annotList[0]=="none":
			notAnnotated+= 1
			continue
		oList.append(o)
		oDict[o] = annotList
		for a in annotList:
			try:
				aDict[a].append(o)
			except KeyError:
				aDict[a] = [o]

print("{} total annotations found accross".format(len(aList)))
print("{} total orthologs".format(totalOrtho))
print("{} had uniprot annotations".format(len(oList)))
print("{} were did not".format(notAnnotated))



#aDict = {uniprotTag : [ortholog1, ortholog2]}


print("\nReading through idmappings...")
#read through the idmappings
#append build dictionary with orthologs keyed to lists of the GO terms
#for all associated annotations
goDict = {}
hadGoList = []
foundAnnotList = []
hadGo = 0
lineCount = 0
with open(mappings, 'r') as infile:
	for line in infile:
		lineCount += 1
		if lineCount % 1e6 == 0:
			print("{} lines processed...".format(lineCount))
			print("\t{} out of {} annotated orthologs found GO terms".format(hadGo, len(oList)))
			print("\t{} out of {} annotations have been found in the mappings file".format(len(foundAnnotList), len(aList)))
			# for x in hadGoList:
			# 	print("{}:\t{}".format(x, ";".join(goDict[x])))
			# exit()
		#noticed that it gets to end after millionth line
		if lineCount == 3e6:
			break
		line=line.strip("\n").split("\t")
		a = line[0]
		try:
			orthos = aDict[a]
			foundAnnotList.append(a)
		except KeyError:
			#This annotation did not have an ortholog associated with it (most)
			continue
		goList = "".join(line[6].split()).split(';')
		if goList[0]=='':
			continue
		for o in orthos:
			try:
				goDict[o] += goList
			except KeyError:
				goDict[o] = goList
				hadGoList.append(o)
				hadGo += 1

print("\n{} Orthologs had GO annotations".format(hadGo))
outName = infileName + "_GO.tsv"
print("\nWriting out results to {}".format(outName))
with open(outName, 'w') as out:
	for o in hadGoList:
		goList = goDict[o]
		uGoList = []
		for go in goList:
			if go not in uGoList:
				uGoList.append(go)
		out.write("{}\t{}\n".format(o, ";".join(uGoList)))



		




