#ortholog_tables
#This directory contains the ortholog results


#contents:
FastOrthoOutput_run4.end
	Initial orthologous groups as output from FastOrtho
singleCopyOrthos.txt
	Initial set of single-copy orhtologs
singleCopyOrthos_revised.txt
	The revised set checked for monophyly of the clades in the checkGroups.txt file (metadata)
	This is considered the best set, and was used for all analyses.
singleCopyAnnotations.tsv
	The annotations file based on the SwissProt and Pfam hits from the Transdecoder pipeline
singleCopyAnnotations.tsv_GO.tsv
	The Gene Ontology annotations for the orhologs, in format for use in GO_MWU.R