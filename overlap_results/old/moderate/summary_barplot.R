

setwd("/Users/grovesdixon/gitreps/convergent_evo_coral/overlap_results/moderate/")
load("objects_for_barplot.Rdata")
#========= plot bargraph similar to Foote ============#
head(tsig)     #branch site significance results
head(asig)     #branch site significance results
head(bs.genes)
head(tt.base)  #the set of convergent substitutions
head(ta.base)  #the set of convergent substitutions between verticals and horizontals
head(sdat)     #the significant sites from branch-site tests
vsdat = sdat[!grepl("Anti", sdat$clade),]
unique(vsdat$clade)
tPairs = sort(unique(paste(tt.base$b1Clade, tt.base$b2Clade, sep="_")))
length(tPairs)



counts = table(tt.base$site)
counts[counts>2]

#gather two-way data
twoWay = counts[counts==1]
dat2 = tt.base[tt.base$site %in% names(twoWay),]
tPairs = sort(unique(paste(tt.base$b1Clade, tt.base$b2Clade, sep="_")))
dat2$pair = paste(dat2$b1Clade, dat2$b2Clade, sep="_")
baseCount = c()
bsCount = c()
siteCount = c()
for (tp in tPairs){
	sub=dat2[dat2$pair==tp,]
	sub$bsSig = sub$ortho %in% bs.genes
	sub$siteSig = sub$site %in% vsdat$site
	baseCount = append(baseCount, length(unique(sub$site)))
	bsCount = append(bsCount, sum(sub$bsSig))
	siteCount = append(siteCount, sum(sub$siteSig))
}
twoRes = data.frame(tPairs, baseCount, bsCount, siteCount)


#gather the three-way results
threeWay = counts[counts==3]
dat3 = tt.base[tt.base$site %in% names(threeWay),]  #all convergent subs that occured accross three species
dat3$bsSig = dat3$ortholog %in% bs.genes
dat3$siteSig = dat3$site %in% vsdat$site
udat3 = dat3[!duplicated(dat3$site),]
a = length(udat3$site) #total three-way substitutions
b = sum(udat3$bsSig)           #total three-way substitutions that were also branch-site significant
c = sum(udat3$siteSig)         #total three-way substitutions that the site was significant in a branch site test
threeRes = c("three-way", a, b, c)


fourWay = counts[counts==4]
length(twoWay) #none


#assemble all together
group = append('three-way', tPairs)
base = append(a, baseCount)
site = append(c, siteCount)
tdf=data.frame(group, base, site)
tdf=tdf[order(tdf$base),]
head(tdf)
t = rbind(tdf$site, tdf$base)
rownames(t) = c('site', 'base')
colnames(t) = tdf$group
col1 = rep('firebrick', nrow(tdf))
col2 = rep('dodgerblue', nrow(tdf))
colors = rbind(col1, col2)
t

#plot absolute values
barplot(t, col=colors, ylab="Molecular Convergence Counts")
legend("topleft", c('convergent', 'convergent + site significant'), fill=c(col1[2], col2[1]))



#normalize by the number of orthologs the groups appear in 
odat = read.table("~/gitreps/convergent_evo_coral/ortholog_tables/singleCopyOrthos_revised.txt", stringsAsFactors=F)
cdat = read.table("~/gitreps/convergent_evo_coral/metadata/vertical_transmitter_clades.txt", stringsAsFactors=F)
colnames(odat) = c('ortho', 'seq')
colnames(cdat) = c('spp', 'clade')
spp=sapply(seq, function(x) strsplit(x, "_")[[1]][1])
odat$spp = spp
head(odat)



#iterate through clades and collect the number of genes they appear in
clades = unique(cdat$clade)
clades = clades[!grepl("Anti", clades)]
nOrthos = c()
clade = clades[1]
pair = tPairs[1]
for (pair in tPairs){
	clade1 = strsplit(pair, "_")[[1]][1]
	clade2 = strsplit(pair, "_")[[1]][2]
	clade1.spps = cdat$spp[cdat$clade==clade1]
	clade2.spps = cdat$spp[cdat$clade==clade2]
	c1sub = odat[odat$spp %in% clade1.spps,]
	c2sub = odat[odat$spp %in% clade2.spps,]
	both = append(c1sub$ortho, c2sub$ortho)
	nOrthos = append(nOrthos, length(unique(both)))
}
rPairs = append("three-way", tPairs)
n.ortho = append(mean(nOrthos), nOrthos)
nres = data.frame(rPairs, n.ortho)
rownames(nres) = nres$rPairs
nres=nres[as.character(tdf$group),]
nres$rPairs==as.character(tdf$group)
tdf$nOrtho = nres$n.ortho
tdf$norm.base = tdf$base / tdf$nOrtho
tdf$norm.site = tdf$site / tdf$nOrtho

t = rbind(tdf$norm.site, tdf$norm.base)
rownames(t) = c('site', 'base')
colnames(t) = tdf$group
col1 = rep('firebrick', nrow(tdf))
col2 = rep('dodgerblue', nrow(tdf))
colors = rbind(col1, col2)
t

#plot absolute values
barplot(t, col=colors, ylab="Molecular Convergence Counts")
legend("topleft", c('convergent', 'convergent + site significant'), fill=c(col1[2], col2[1]))








