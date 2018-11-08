
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.


# Edit these to match your data file names: 
setwd("~/gitreps/convergent_evo_coral/overlap_results/moderate/GO_MWU")

# #select input file
# input="GOMWU_input1_branchSitesVertical.csv"
# input="GOMWU_input2_convergence_overlap.csv"
# input="GOMWU_input3_bs_convergence_overlap.csv"
# input="GOMWU_input3_bs_convergence_overlapHH.csv"


# goAnnotations="singleCopyAnnotations.tsv_GO.tsv"   #old version made from idmapping_selected.tab (don't use this)
goAnnotations="singleCopyAnnotations_GO_gomwu.tsv" #new better version made from goa_uniprot_all.gaf
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
# goDivision="CC"     # either MF, or BP, or CC
source("gomwu.functions.R")


#SET UP SETS OF VARS FOR LOOPED RUNS
ABS.VALUE=0.99

inputs = c(
	'GOMWU_input1_branchSitesVertical.csv',
	'GOMWU_input2_covergenceVV.csv',
	'GOMWU_input3_bs_convergence_overlapVV.csv',
	'GOMWU_input1_branchSitesHorizontal.csv',
	'GOMWU_input2_covergenceHH.csv',
	'GOMWU_input3_bs_convergence_overlapHH.csv',
	'GOMWU_input4_flagged_convergence_overlapVV.csv',
	'GOMWU_input4_flagged_convergence_overlapHH.csv'
	)


goDivisions = c('CC', 'MF', 'BP')
goDivisions=c('CC')


#Loop through inputs and GO divisions running GO_MWU for all
for (input in inputs){
	for (goDivision in goDivisions){
		print('----------------')
		print(paste(paste('RUNNING FOR ', input), goDivision))
		print('----------------')

		#Run GO MWU
		gomwuStats(input, goDatabase, goAnnotations, goDivision,
			perlPath="perl", 
			largest=0.1,  
			smallest=200,   
			clusterCutHeight=0.25, 
			Alternative="g"
		)
		
		#Save top 20
		resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
		res=read.table(resName, header = T)
		res=res[order(res$pval),]
		tabOut = paste('..', sub(".csv", "_top20.tsv", resName), sep="/")
		tabOut
		write.table(head(res, n=20), file=tabOut, quote=F, sep="\t")
	}
}


# Plotting results
quartz()
gomwuPlot(input,goAnnotations,goDivision,
	absValue=ABS.VALUE,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
	level2=0.05, # FDR cutoff to print in regular (not italic) font.
	level3=0.01, # FDR cutoff to print in large bold font.
	txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5, # height of the hierarchical clustering tree
#	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

#Use GO_MWU_viewResults.R to look at genes associated with these

