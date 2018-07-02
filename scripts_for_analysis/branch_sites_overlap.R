#branch_sites_overlap.R
#Groves Dixon
#last updated 7-1-18

setwd("~/gitreps/convergent_evo_coral")
source('scripts_for_analysis/coral_reproduction_functions.R')


#---- CUTOFF SELECTION ----#
#The sets of cutoffs are stored in separate R scripts

source("scripts_for_analysis/moderate_cutoff_set.R") 

# source("scripts_for_analysis/relaxed_cutoff_set.R")

#------ LOAD DATA ------#

#read in lrt results from each the branch sites tests
#target datasets
gdat = read.table("branch_sites_tests/bs_lrt_Gacrhelia.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)
mdat = read.table("branch_sites_tests/bs_lrt_Montipora.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)
pdat = read.table("branch_sites_tests/bs_lrt_Porites.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)
qdat = read.table("branch_sites_tests/bs_lrt_Pocillopora.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)
vdat = read.table("branch_sites_tests/bs_lrt_allVertical.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)



#"anti" target datasets
a.gdat = read.table("branch_sites_tests/bs_lrt_anti_Gacrhelia.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)
a.mdat = read.table("branch_sites_tests/bs_lrt_anti_Montipora.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)
a.pdat = read.table("branch_sites_tests/bs_lrt_anti_Porites.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)
a.qdat = read.table("branch_sites_tests/bs_lrt_anti_Pocillopora.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)
a.vdat = read.table("branch_sites_tests/bs_lrt_anti_allVertical.tsv", sep = "\t", header = T, quote="", stringsAsFactors=F)

#---- FORMAT AND MERGE ----#

#specify columns
colnames(gdat) = paste('g', colnames(gdat), sep = "_")
colnames(mdat) = paste('m', colnames(mdat), sep = "_")
colnames(pdat) = paste('p', colnames(pdat), sep = "_")
colnames(qdat) = paste('q', colnames(qdat), sep = "_")
colnames(vdat) = paste('v', colnames(vdat), sep = "_")


colnames(a.gdat) = paste('g', colnames(a.gdat), sep = "_")
colnames(a.mdat) = paste('m', colnames(a.mdat), sep = "_")
colnames(a.pdat) = paste('p', colnames(a.pdat), sep = "_")
colnames(a.qdat) = paste('q', colnames(a.qdat), sep = "_")
colnames(a.vdat) = paste('v', colnames(a.vdat), sep = "_")

### add merging column
#targets
gdat$prot=gdat$g_contig
mdat$prot=mdat$m_contig
pdat$prot=pdat$p_contig
qdat$prot=qdat$q_contig
vdat$prot=vdat$v_contig

#antis
a.gdat$prot=a.gdat$g_contig
a.mdat$prot=a.mdat$m_contig
a.pdat$prot=a.pdat$p_contig
a.qdat$prot=a.qdat$q_contig
a.vdat$prot=a.vdat$v_contig

#merge together
#targets
c0=merge(gdat, mdat, by = 'prot', all=T)
c1=merge(c0, pdat, by='prot', all=T)
c2 = merge(c1, qdat, by = 'prot', all=T)
cdat = merge(c2, vdat, by = 'prot', all=T)
rownames(cdat) = cdat$prot
head(cdat)

#antis
a.c0=merge(a.gdat, a.mdat, by = 'prot', all=T)
a.c1=merge(a.c0, a.pdat, by='prot', all=T)
a.c2 = merge(a.c1, a.qdat, by = 'prot', all=T)
a.cdat = merge(a.c2, a.vdat, by = 'prot', all=T)
rownames(a.cdat) = a.cdat$prot
head(a.cdat)

#total numbers of genes with bs data
nrow(cdat)
nrow(a.cdat)


#----- SIGNIFICANCE OVERLAP -----#

#look for genes with significance in the branch sites results
print(paste("Number of significant lineages must be at least", bs.min.lineage))
bs.genes = bs_sig_overlap(df=cdat, ptype=bs.ptype)      #significant in at least one vertical transmitter test
a.bs.genes = bs_sig_overlap(df=a.cdat, ptype=bs.ptype)  #signifiant in at least one horizontal tramsitter test
all.bs.genes = unique(c(bs.genes, a.bs.genes))          #combined


#build matrix with booleans for whether each test was significant
tsig = make_sig_matrix(cdat, bs.ptype, bs.alpha)
asig = make_sig_matrix(a.cdat, bs.ptype, bs.alpha)
colnames(tsig) = paste('target', colnames(tsig), sep="_")
colnames(asig) = paste('anti', colnames(asig), sep="_")
sig =merge(tsig, asig, by = 0)
rownames(sig) = sig$Row.names
sig$Row.names<-NULL
head(sig) 


#OUTPUT GO ENRICHMENT FILE
#set up binary for whether each gene was significant in at least one test for vertical transmitters
genes = unique(c(rownames(tsig), c(rownames(asig))))
score = rep(0, length(genes))
score[genes %in% bs.genes]<-1
goin1 = data.frame(genes, score)
colnames(goin1) = c("gene", "score")
write.csv(goin1, file=paste(results.path, "GO_MWU/GOMWU_input1_branchSites.csv",sep='/'), quote=F, row.names=F)

#breakdown of scores:
table(score)
plot(density(score, na.rm=T))


#--- ADD CONVERGENCE DATA ---#

#read in the ancestral convergence results
ll = load(ancestral.file) #This is the ancestral results file generated with parse_multi_rstV2.py.
head(dat)
dim(dat)

#apply the confidence cutoff to the ancestral states
dat = dat[dat$b1post > anc.post.cutoff & dat$b2post > anc.post.cutoff,]
dim(dat)
dat$site = paste(dat$ortholog, dat$b1pos, sep="_")
anc.genes = unique(dat$ortholog)


#===== OUTPUT GO ENRICHMENT FOR CONVERGENCE ALONE =====#

#get counts for target-target (tt) overlap and 
#subtract those for target-antiTarget (ta)
#the 'score' in the output table is the difference in the 
#number of convergent substitutions that occured between
#pairs from the target group, and pairs from the nontarget group

#logic is that the ta group demonstrates a sort of null,
#which we compare the frequency of convergence events in
#the tt group against. Only GO terms that have convergence
#preferentially in the tt group will be enriched.


#separate the convergence events based on traits
c1 = grepl("Anti", dat$b1Clade)
c2 = grepl("Anti", dat$b2Clade)
conv = dat$genConvergent
tt.base = dat[!c1 & !c2 & conv, ]
aa.base = dat[c1 & c2 & conv, ]
ta1 = dat[!c1 & c2 & conv,]
ta2 = dat[c1 & !c2 & conv,]
ta.base = rbind(ta1, ta2)

#sum up the number of convergence events for each gene
tt.sums = data.frame(table(tt.base$ortholog))
ta.sums = data.frame(table(ta.base$ortholog))
aa.sums = data.frame(table(aa.base$ortholog))

#merge and output for GO enrichment
m = merge(tt.sums, ta.sums, by = 'Var1')
gene = as.character(m$Var1)
score = m$Freq.x
score[score>1]<-1
rest.of.genes = genes[!genes %in% gene]
rest.of.scores = rep(0, length(rest.of.genes))
go.genes = append(gene, rest.of.genes)
go.scores = append(score, rest.of.scores)
goin2 = data.frame(go.genes, go.scores)
colnames(goin2) = c('gene', 'score')
write.csv(goin2, file=paste(results.path, "GO_MWU/GOMWU_input2_convergence_overlap.csv",sep='/'), quote=F, row.names=F)
plot(density(goin2$score));abline(v=0)
sum(goin2$score)
#=========================================

#SUBSET FOR THE MINIMUM NUMBER OF CONVERTENT SITES REQUIRED

counts = table(dat$site)
multiSites = names(counts[counts>=min.convergent])
length(multiSites)
subMulti = dat[dat$site %in% multiSites,]
nrow(subMulti)
length(unique(subMulti$ortholog))


#SUBSET FOR CONVERGENT AND SIGNIFICANT FOR BRANCH SITES
#These are the convergent substitutions that occured in 
#genes that were significant in the branch sites tests

sub2 = subMulti[subMulti$ortholog %in% all.bs.genes,]
sub2=sub2[order(sub2$site),]
print(paste("Total overlapping substitutions among genes that were in the branch sites overlap set =", nrow(sub2)))
nrow(sub2)
print("Total genes:")
length(unique(sub2$ortholog))


bs.conv.subs = sub2[sub2$genConvergent,]
print(paste("Total convergent substitutions in overlapset =", nrow(bs.conv.subs)))
print(paste("Total genes with convergent subtitutions =", length(unique(bs.conv.subs$ortholog))))


#SPECIFY FOR CONVERGENCE AND BRANCH SITES SIGNIFICANCE WITH REGARD TO TRAIT OF INTEREST

#tt = target-target (ie vertical-vertical)
##ta = target-antiTarget
##aa = anti-anti

t1 = !grepl('Anti', bs.conv.subs$b1Clade)
t2 = !grepl('Anti', bs.conv.subs$b2Clade)

tt = bs.conv.subs[t1 & t2,]
ta = bs.conv.subs[(t1 & !t2) | (!t1 & t2),]
aa = bs.conv.subs[!t1 & !t2,]

#remove sister comparisons from ta
sister1 = sub("Anti", "", ta$b1Clade)
sister2 = sub("Anti", "", ta$b2Clade)


#print out counts for the overlapping branch sites and convergence results

print("Numbers of substitutions from each trait combination:")
nrow(tt)
nrow(ta)
nrow(aa)

print("Numbers of unique genes for each trait combination:")
length(unique(tt$ortholog))
length(unique(ta$ortholog))
length(unique(aa$ortholog))


print("All subs accounted for?")
sum(c(nrow(tt), nrow(ta), nrow(aa)))==nrow(bs.conv.subs)

#===== OUTPUT GO ENRICHMENT FOR BS-CONVERGENCE OVERLAP =====#

#output for GO enrichment
head(goin1)
head(goin2)
conv.genes = goin2$gene[goin2$score==1]
print(paste("Total genes with at least one convergence event =", length(conv.genes)))
print(paste("Total genes with at least one significant branch site test =", length(bs.genes)))
genes = goin1$gene[goin1$gene %in% goin2$gene]
score = rep(0, length(genes))
score[genes %in% bs.genes & genes %in% conv.genes]<-1
goin3 = data.frame(genes, score)
colnames(goin3) = c('gene', 'score')
write.csv(goin3, file=paste(results.path, "GO_MWU/GOMWU_input3_bs_convergence_overlap.csv",sep='/'), quote=F, row.names=F)
plot(density(goin3$score)); abline(v=0)
sum(goin3$score)
#====================




#--------- SUBSET FOR SITES FLAGGED IN THE BRANCH-SITES OUTPUTS ------------#
#build these files with extract_sig_bs_positions.py
#targets
sdat1 = read.table("branch_sites_tests/gachrelia_sigBsSites.tsv", header = T)
sdat2 = read.table("branch_sites_tests/montipora_sigBsSites.tsv", header = T)
sdat3 = read.table("branch_sites_tests/porites_sigBsSites.tsv", header = T)
sdat4 = read.table("branch_sites_tests/pocillopora_sigBsSites.tsv", header = T)
sdat5 = read.table("branch_sites_tests/allVertical_sigBsSites.tsv", header = T)

#antis
sdat6 = read.table("branch_sites_tests/anti_gachrelia_sigBsSites.tsv", header = T)
sdat7 = read.table("branch_sites_tests/anti_montipora_sigBsSites.tsv", header = T)
sdat8 = read.table("branch_sites_tests/anti_porites_sigBsSites.tsv", header = T)
sdat9 = read.table("branch_sites_tests/anti_pocillopora_sigBsSites.tsv", header = T)
sdat10 = read.table("branch_sites_tests/anti_allVertical_sigBsSites.tsv", header = T)

#set the clade names for the datasets
#targets
sdat1$clade='Galaxia'
sdat2$clade="Montipora"
sdat3$clade='Porites'
sdat4$clade='Pocilloporid'
sdat5$clade="AllVertical"

#antis
sdat6$clade='AntiGalaxia'
sdat7$clade="AntiMontipora"
sdat8$clade='AntiPorites'
sdat9$clade='AntiPocilloporid'
sdat10$clade='AntiAllVertical'


#assmble into single dataframe
sdat = rbind(sdat1, sdat2, sdat3, sdat4, sdat5, sdat6, sdat7, sdat8, sdat9, sdat10)
sdat$ortholog = sub("_ALT.codeml", "", sdat$file)
sdat$site = paste(sdat$ortholog, sdat$position, sep="_")
sdat=sdat[sdat$posterior>=bs.site.post.cutoff,]
head(sdat)


#subset for site overlap between flagged sites and convergent subs among the branch-sites significant genes
osub = bs.conv.subs[bs.conv.subs$site %in% sdat$site,]
osite = sdat[sdat$site %in% bs.conv.subs$site,]
head(osub)
head(osite)
print(paste("Total convergent sites also flagged for positive selection =", nrow(osub)))
print(paste("Total unique genes =", length(unique(osub$ortholog))))


#SUBSET FOR CONVERGENT SUBS FLAGGED IN SAME CLADE
#This includes the 'allVertical' and 'allAntiVertical'

#merge the two datasets together
m = merge(osub, sdat, by = 'site')
u=unique(c(m$b1Clade, m$b2Clade))
u2 = unique(m$clade)
sum(u2 %in% u)==length(u)

#set up the convergent subs flagged in same clade
m1 = m[m$b1Clade==m$clade | m$b2Clade==m$clade,]

#set up convergent subs flagged in the appropritate 'all' tests
con.b1Target = !grepl("Anti", m$b1Clade) #the convernt branch 1 is target
con.b2Target = !grepl("Anti", m$b2Clade) #the convergent branch 2 is target
flagClade = m$clade

#gather convergent sites also flagged in the 'all' branch sites tests
all.target = m[con.b1Target & con.b2Target & flagClade=='AllVertical',]
all.ta =  m[con.b1Target & !con.b2Target & flagClade=='AllVertical',]
all.at =  m[!con.b1Target & con.b2Target & flagClade=='AllVertical',]
all.anti = m[!con.b1Target & !con.b2Target & flagClade=='AntiAllVertical',]

#assemble together
m2 = rbind(m1, all.target, all.anti, all.ta, all.at)

#print summary of results
nrow(all.target)
nrow(all.ta)
nrow(all.at)
nrow(all.anti)

head(m2)
nrow(m2)
print(paste("Total convergent subs flagged for positive selection in one of the convergent lineages=", nrow(m2)))
print(paste("Total unique genes =", length(unique(m2$ortholog.x))))


#now subset for tt, ta, or aa convergence pairs
m2$t1 = !grepl('Anti', m2$b1Clade)
m2$t2 = !grepl('Anti', m2$b2Clade)
head(m2)
tt= m2[m2$t1 & m2$t2,]
ta = m2[ (m2$t1 & !m2$t2) | (!m2$t1 & m2$t2),]
aa = m2[!m2$t1 & !m2$t2,]

print("Total final convergent substitutions that meet all criteria:")
print(paste("TT =", length(tt$ortholog.x)))
print(paste("TA =", length(ta$ortholog.x)))
print(paste("AA =", length(aa$ortholog.x)))

print("Total final convergent genes that meet all criteria:")
print(paste("TT =", length(unique(tt$ortholog.x))))
print(paste("TA =", length(unique(ta$ortholog.x))))
print(paste("AA =", length(unique(aa$ortholog.x))))


#how many convergence events occured accross more than one pair of clades?
ttc = table(table(tt$site))
tac = table(table(ta$site))
table(tt$site)[table(tt$site)>1]


#---- FORMAT FINAL RESULTS AND OUTPUT -----#

#FINAL SET OF SUBTITUTIONS FITTING ALL CRITERIA
annots = read.table("ortholog_tables/singleCopyAnnotationsWithDescriptions.tsv", header = T, sep="\t")
tt.out = paste(results.path, 'verticalConvergentRes.tsv', sep="/")
annots$ortholog.x = annots$ortho
ttm = merge(tt, annots, by = "ortholog.x", all.x=T)
ttm=ttm[!colnames(ttm) %in% c('conv', 'divergent', 'para', 'none', 'genConvergent', 'genNone', 'file', 'ortholog.y', 't1', 't2')]
head(ttm)
nrow(ttm)
write.table(ttm, file= tt.out, sep="\t", quote=F, row.names=F)


#final unique genes (use these to output trees and intermediate files)
genes.out = paste(results.path, 'verticalConvergentGenes.txt', sep="/")
write.table(unique(tt$ortholog.x), file= genes.out, sep="\t", quote=F, row.names=F, col.names=F)


#============ FINAL GO OUTPUT ===========#
#This one simply outputs the total number of convergent substitutions also 
#flagged as positively selected in the appropriate lineage (including allVertical)

genes = unique(c(rownames(tsig), c(rownames(asig))))
fscores = data.frame(table(ttm$ortholog.x))
colnames(fscores) = c('gene', 'score')
rest.of.genes = genes[!genes %in% fscores$gene]
rest.of.scores = rep(0, length(rest.of.genes))
rest.df = data.frame(rest.of.genes, rest.of.scores)
colnames(rest.df) = c('gene', 'score')
finalGo = rbind(fscores, rest.df)
finalGo$score[finalGo$score>1]<-1
head(finalGo)
dim(finalGo)
table(finalGo$score)
goout.path = paste(results.path, "GO_MWU/GOMWU_input4_final.csv", sep="/")
write.csv(finalGo, file= goout.path, quote=F, row.names=F)


#save objects for building barplot
save(bs.genes, tsig, asig, tt.base, ta.base, sdat, file=paste(results.path, "objects_for_barplot.Rdata", sep="/"))


#------------------------------------------- old stuff below -------------------------------------------#









#subset for clades
fsub = osub[osub$b1Clade %in% clades & osub$b2Clade %in% clades,]

#format finalized set
osite = osite[,c('site', 'clade', 'posterior')]
colnames(osite) = c('site', 'siteStat.clade', 'sitePosterior')
fm = merge(osite, fsub, by = 'site')
annots$ortholog=rownames(annots)
fm2 = merge(fm, annots, by = 'ortholog')
fm2
dim(fm2)
write.table(fm2, file=paste(results.path, 'final_overlap_results.tsv', sep='/'), quote=F, row.names=F, sep="\t")

#--- for gene IDs
sp = fm2$swissProt
o = fm2$ortholog
spString = sp[1]
ortho=o[1]
genes = unlist(strsplit(spString, ";"))
sps = sapply(genes, function(x) strsplit(x, "|", fixed=T)[[1]][1])
spres = data.frame(rep(ortho, n=length(sps)), sps)
colnames(spres) = c('ortholog', 'sp')

for (i in 2:length(sp)){
	spString = sp[i]
	ortho=o[i]
	genes = unlist(strsplit(spString, ";"))
	sps = sapply(genes, function(x) strsplit(x, "|", fixed=T)[[1]][1])
	sub = data.frame(rep(ortho, n=length(sps)), sps)
	colnames(sub) = c('ortholog', 'sp')
	sub=sub[sub$sp !='none',]
	spres=rbind(spres, sub)
	rownames(spres) = 1:nrow(spres)
}
spres
write.table(spres, file=paste(results.path, 'swissProt_table.tsv', sep="/"), quote=F, row.names=F, sep="\t")
write.table(spres$sp, file=paste(results.path, 'swissProt_ids.tsv', sep="/"), quote=F, row.names=F, sep="\t", col.names=F)


#to gather gene trees and amio acid alignments see section OUTPUT FILES FOR SELECTED GENES in transdecorder_fastortho.txt













#output for GO
all.genes = rownames(cdat)
score = rep(0, length(all.genes))
score[all.genes %in% overconv$geneId]<-1
goin = data.frame(gene=all.genes, score=score)
head(goin)
write.csv(goin, file="GO_MWU/overlap_counts.csv", quote=F, row.names=F)






bs=length(bs.genes)
not.bs = nrow(cdat) - bs
cv = length(conv.genes)
not.cv = length(anc.genes) - cv

tab = rbind(c(bs, not.bs), c(cv, not.cv))
colnames(tab) = c('in', 'out')
rownames(tab) = c('bs', 'conv')












d=mdat
p.type = 'adj.p'
p.type = 'p.values'
CUT=0.1

get_overlap = function(d, condat, p.type, CUT, FOREGROUND){
	s=d[d[,p.type]<CUT,]
	s$geneId = s$contig
	print(paste(nrow(s), 'significant genes'))
	o=condat[condat$geneId %in% s$contig,]
	print(paste('overlapping sites =', nrow(o)))
	m=merge(o, s, by = 'geneId')
	print(paste(length(unique(m$protein.name)), 'genes:'))
	for(g in unique(m$protein.name)){
		print(g)
	}
	if (nrow(m) > 0){
		m$foreground=FOREGROUND
		return(m)
	}
}
	
g.o = get_overlap(gdat, condat, p.type, CUT, 'Gacrhelia')
m.o = get_overlap(mdat, condat, p.type, CUT, 'Montipora')
p.o = get_overlap(pdat, condat, p.type, CUT, 'Porites')
v.o = get_overlap(vdat, condat, p.type, CUT, 'vertical')
all = rbind(g.o, m.o, p.o, v.o)
all=all[order(all$p.values),]
# all=all[!duplicated(all$gene_pos),]
u.all = all[!duplicated(all$protein.name),]
u.all
nrow(u.all)
write.table(u.all, "results/selected_genes/unique_convergent_bs_genes.tsv", sep = "\t", row.names = F, quote=F)
write.table(all, "results/selected_genes/all_convergent_bs_genes.tsv", sep = "\t", row.names = F, quote=F)

#export the data for GO and KOGG enrichment tests using MWU-tests
out2 = data.frame(out$contig, -log(out$p.values, 10))
colnames(out2) = c('gene', 'logp')
head(out2)
nrow(out2)
#write out for GO
#use this file as input for 
write.table(out2, 'branch_site_LRT_results_for_GOmwu.csv', quote = F, row.names = F, sep = ",")
ps = out$p.values
x = ps < 0.05
ps[x==TRUE]<-1
ps[x==FALSE]<-0
out3 = data.frame(out$contig, ps)
colnames(out3) = c('isogroup', 'sig')
write.table(out3, 'branch_site_LRT_results_for_Fisher_GOmwu.csv', quote = F, row.names = F, sep = ",")

#write out for KOGG
write.table(out2, 'branch_site_LRT_results_for_KOGGmwu.csv', quote = F, row.names = F, sep = ",")









#---------- build the figure 2 barplot from Foote ----------#
#groupings are: combined vertical transmitters, Gacrhelia and Montipora, Gacrhelia and Porites, Montipora and Porites
all.bs = rbind(gdat, mdat, pdat, vdat)
all.sig = all.bs[all.bs$p.values < 0.1,]
all.sig = all.sig[order(all.sig$p.values),]
all.sig = all.sig[!duplicated(all.sig$contig),]
sig.genes = all.sig$contig
head(all.sig)
nrow(all.sig)


#get combined vertical transmitter gene counts
head(condat)
uni.condat = condat[!duplicated(condat$geneId),]
tot.genes = nrow(uni.condat)
bs.genes = uni.condat[uni.condat$geneId %in% sig.genes,]
bs.genes
tot.bs = nrow(bs.genes)
tot.non.bs = tot.genes - tot.bs
combined = c(tot.bs, tot.non.bs)

#get combined vertical transmitter site counts
uni.condat.sites = condat[!duplicated(condat$gene_pos),]
tot.sites = nrow(uni.condat.sites)
bs.sites = uni.condat.sites[uni.condat.sites$geneId %in% sig.genes,]
tot.bs.sites = nrow(bs.sites)
tot.non.sites = tot.sites-tot.bs.sites
combined.sites=c(tot.bs.sites, tot.non.sites)

#lazily done funciton to get pair-wise data for barplot
#!NOTE this only works if the different clades compared have unique letters
get_pair_data = function(abrv1, abrv2, abrv3="None", sites=FALSE){
	#subset for parallel aa substitutions
	sub1 = dat[grep(abrv1, dat$paraPairs),]
	sub1 = sub1[grep(abrv2, sub1 $paraPairs),]
	#subset for convergent aa substitutions
	sub2 = dat[grep(abrv1, dat$convPairs),]
	sub2 = sub2[grep(abrv2, sub2 $convPairs),]
	npara = nrow(sub1)
	nconv = nrow(sub2)
	comb0 = rbind(sub1, sub2)
	uni.sites = comb0[!duplicated(comb0$gene_pos),]
	tot.sites = nrow(uni.sites)
	uni.genes = comb0[!duplicated(comb0$geneId),]
	tot.genes = nrow(uni.genes)
	bs.comb = uni.genes[uni.genes$geneId %in% sig.genes,]
	tot.bs = nrow(bs.comb)
	tot.non.bs = tot.genes - tot.bs
	bs.sites = uni.sites[uni.sites$geneId %in% sig.genes,]
	tot.bs.sites = nrow(bs.sites)
	print(paste(npara, 'parallel substitutions found'))
	print(paste(nconv, 'convergent substitutions found'))
	print(paste(tot.sites, 'total sites found with molecular convergence'))
	print(paste(tot.genes, 'total genes with molecular convergence'))
	print(paste(tot.genes, 'total genes with molecular convergence'))
	print(paste(tot.bs, 'total genes with molecular convergence showed evidence of positive selection in at least one branch-sites test'))
	if (sites){
		print("Returning Site Count Data")
		return(c(tot.bs.sites, tot.sites-tot.bs.sites))
	}
	else{
		print("Returning Gene Count Data")
		return(c(tot.bs, tot.non.bs))
	}
	
}

##### Build gene-wise barplot

#get pairwise convergence data
g.m = get_pair_data(abrv1 = "G", abrv2 = "M")
g.p = get_pair_data(abrv1 = "G", abrv2 = "P")
m.p = get_pair_data(abrv1 = "M", abrv2 = "P")


#assemble data for barplot
res=rbind(combined, g.m, g.p, m.p)
colnames(res) = c('bs', 'conv')
res=res[order(res[,2]),]
bp.dat = t(res)
labs = c('Combined', 'G.acrhelia and Montipora', 'G.acrhelia and Porites', 'Montipora and Porites')
labs = c('Combined', 'G.acrhelia and Montipora', 'G.acrhelia and Porites', 'Montipora and Porites')

gene.res = res

#plot
barplot(bp.dat, col=c('firebrick', 'deepskyblue3'), space=0, names = labs, ylab="Number of Genes Showing Molecular Convergence")
save(bp.dat, labs, file='figures/genes_count_barplot_data.Rdata')


##### Build site-wise barplot
#get pairwise convergence data
g.m = get_pair_data(abrv1 = "G", abrv2 = "M", sites=TRUE)
g.p = get_pair_data(abrv1 = "G", abrv2 = "P", sites=TRUE)
m.p = get_pair_data(abrv1 = "M", abrv2 = "P", sites=TRUE)


#assemble data for barplot
res=rbind(combined.sites, g.m, g.p, m.p)
colnames(res) = c('bs', 'conv')
res=res[order(res[,2]),]
bp.dat = t(res)
labs = c('Combined', 'G.acrhelia and Montipora', 'G.acrhelia and Porites', 'Montipora and Porites')
labs = c('Combined', 'G.acrhelia and Montipora', 'G.acrhelia and Porites', 'Montipora and Porites')

#plot
barplot(bp.dat, col=c('firebrick', 'deepskyblue3'), space=0, names = labs, ylab="Number of Sites Showing Molecular Convergence")
save(bp.dat, labs, file='figures/sites_count_barplot_data.Rdata')

#plot site by site
res
gene.res


