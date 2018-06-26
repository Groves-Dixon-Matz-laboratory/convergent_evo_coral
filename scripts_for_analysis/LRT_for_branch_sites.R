#LRT_for_branch_sites_models
#Groves Dixon
#10/20/17
#R script LRT_for_branch_site_models.R


#set working directory
rm(list=ls())
setwd("~/gitreps/coral_reproductive_evolution")


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
null.name = paste(paste("results/branch_sites/finalVertical/nullLikelihoods_branchSites", foreground, sep="_"), "tsv", sep=".")
alt.name = paste(paste("results/branch_sites/finalVertical/altLikelihoods_branchSites", foreground, sep="_"), "tsv", sep=".")
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

hist(dat$p.values)



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



#--- old stuff below ---#


#export the data for GO and KOGG enrichment tests using MWU-tests
out2 = data.frame(rdat$contig, -log(rdat$p.values, 10))
colnames(out2) = c('gene', 'logp')
head(out2)
nrow(out2)
#write out for GO
o2.name = paste(paste("gomwu", foreground, sep="/"), "bs_lrt_goInput.csv", sep="_")
write.table(out2, o2.name, quote = F, row.names = F, sep = ",")
#write out for Fisher GO
ps = rdat$p.values
x = ps < 0.05
ps[x==TRUE]<-1
ps[x==FALSE]<-0
out3 = data.frame(rdat$contig, ps)
colnames(out3) = c('isogroup', 'sig')
o3.name = paste(paste("gomwu", foreground, sep="/"), "bs_lrt_Fisher_goInput.csv", sep="_")
write.table(out3, o3.name, quote = F, row.names = F, sep = ",")

#write out for KOGG
write.table(out2, 'branch_site_LRT_results_for_KOGGmwu.csv', quote = F, row.names = F, sep = ",")
