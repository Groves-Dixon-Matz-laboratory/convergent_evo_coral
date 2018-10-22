#LRT_for_branch_sites_models
#Groves Dixon
#10/20/17
#R script LRT_for_branch_site_models.R


#set working directory
rm(list=ls())
setwd("~/gitreps/convergent_evo_coral/branch_sites_tests")


#likelihood ratio test function
lrt = function(la, lo, df){
  G = 2*(la - lo)
  # print(G)
  p.values = pchisq(G, df, ncp = 0, lower.tail = F)
  return(p.values)
}

#choose the branch sites comparison you want:
foreground = "gacrhelia"
foreground = "montipora"
foreground = "porites"
foreground = "pocillopora"
foreground = "allVertical"
foreground = "anti_gacrhelia"
foreground = "anti_montipora"
foreground = "anti_porites"
foreground = "anti_pocillopora"
foreground = "anti_allVertical"




#read in the data from the null and alternative models for the branch-site test for positive selection
#see PAML manual for descriptions of these models
null.name = paste(paste("nullLikelihoods_branchSites", foreground, sep="_"), "tsv", sep=".")
alt.name = paste(paste("altLikelihoods_branchSites", foreground, sep="_"), "tsv", sep=".")
null = read.table(null.name, header = T) #file generated using parse_codeml_branch_sites.py (see Positive_Selection_Walkthrough.txt)
alt = read.table(alt.name, header = T)   #file generated using parse_codeml_branch_sites.py (see Positive_Selection_Walkthrough.txt)


#these should show the contig, the number of parameters for the model, and it's log likelihood
head(null); head(alt)


#run lrt between the models
contig = alt$contig
contigs.null = null$contig
sum(contig == contigs.null) == nrow(alt) #double-check the files were assembled correctly
la = alt$likelihood
lo = null$likelihood
p.values = lrt(la, lo, df=1) #note degrees of freedom is difference in number of parameters between models (in this case always 1)

#record the results of the test in a new dataframe
dat = data.frame(la, lo, p.values, contig, contigs.null)
head(dat)

#doublecheck the contigs match
sum(dat$contig == dat$contigs.null) == nrow(dat)


#get adjusted p values for multiple tests
dat$adj.p = p.adjust(dat$p.values, method = 'BH')

#get an idea of how many genes are significant
cuts = c(0.1, 0.05, 0.01, 0.001, 0.0001)
for (i in cuts){
	sub = dat[dat$adj.p < i,]
	unadjust = dat[dat$p.values < i,]
	print(paste(paste(paste(i, nrow(sub)), nrow(unadjust)), (nrow(sub)/nrow(dat))*100))
}

hist(dat$p.values, main="P-values before final QC")


#WRITE OUT THE RESULTS
#get gene names
gdat = read.table("ortholog_tables/singleCopyAnnotations.tsv", header = T)[,c(1,2)]
colnames(gdat) = c('contig', 'swissProt_hits')
m=merge(dat, gdat, by = 'contig')
m$contigs.null<-NULL
print("All genes accounted for:")
nrow(m) == nrow(dat)
nrow(m) == nrow(alt)
nrow(m)
rdat = m[order(m$p.values),]
nrow(rdat)


#output
out.name = paste(paste("results/branch_sites/finalVertical/bs_lrt", foreground, sep="_"), "tsv", sep=".")
write.table(rdat, file=out.name, row.names = F, quote=F, sep = "\t")


#----- DOUBLE-CHECK THAT RESULTS MATCH WITH DESCRIPTION IN PAML MANUAL:

#from page 31 on the PAML manual:
#"Suppose your 2 delta likelihood = 2.71, you will get 0.10 from chi2, then the p value for the mixture is 0.1/2 = 5%
#We recommend that you use X2 (with critical values 3.84 and 5.99) instead of the mixture to guide against violations of model assumptions"
#the second part, (with critical values 3.84 and 5.99), amounts to 
#not dividing the X2 by 2 to get the pvalue.

#here is an example of running the function above
#showing that 
G = 2.71
p.values = pchisq(G, 1, ncp = 0, lower.tail = F)
p.values #gives 0.1

G = 3.84
p.values = pchisq(G, 1, ncp = 0, lower.tail = F)
p.values #gives 0.05

#so by not dividing by 2 we are following the recommendation





