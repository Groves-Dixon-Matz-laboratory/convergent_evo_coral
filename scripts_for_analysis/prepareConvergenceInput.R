#prepareConvergenceInput.R
#This takes the large output file from parse_multi_rstV2.py
#It divides the file up into manageble pieces for different overlap sets
setwd("~/gitreps/coral_reproductive_evolution/")

#--- UPLOAD CONVERGENCE DATA ---#
#read in the ancestral convergence results
dat = read.table("results/ancestralConvergence/verticalAncestralRes.tsv", sep = "\t", header = T, stringsAsFactors=F, fill=T) #file generated using parse_multi_rstV2.py
dat = read.table("results/ancestralConvergence/ballancedVerticalAncestralRes.tsv", sep = "\t", header = T, stringsAsFactors=F, fill=T) #file generated using parse_multi_rstV2.py


head(dat)
anc.genes = unique(dat$ortholog)



#--- CALL SUBSTITUTION TYPES ---#
#parallel, convergent, divergent, or none
a1=dat$b1Anc
a2=dat$b2Anc
n1=dat$b1New
n2=dat$b2New
same.a = a1==a2
same.n = n1==n2

parallel = same.a & same.n
divergent = same.a & !same.n
convergent = !same.a & same.n
none = !same.a & !same.n

nd = sum(divergent)
nc = sum(convergent)
np = sum(parallel)
nn = sum(none)
x=c(nc, nd, nn, np)
#all types accounted for?
sum(x) == nrow(dat)
type = rep("divergent", nrow(dat))
type[parallel]<-'parallel'
type[convergent]<-'convergent'
type[none]<-'none'
table(type)


#output subdatasets
dat$sub.type = type
dat$conv = type == 'convergent'
dat$para = type == 'parallel'
dat$divergent = type == 'divergent'
dat$none = type == 'none'
dat$genConvergent = dat$conv | dat$para
dat$genNone = dat$divergent | dat$none
head(dat)
dat.backup = dat


#------- SUBSET FOR OVERLAPPING SUBSTITUTIONS WITHIN FULL CLADE LINEAGES ----#
#These have "All" in the branch names
#subset datat for only 'All' lineage subs
subAll = dat[grepl("All", dat$b1Name) & grepl("All", dat$b2Name),]
subAll = subAll[subAll$b1post > 0.8 & subAll$b2post > 0.8, ]
dat = subAll
# save(dat, file="datasets/ancestral_recon_full_clades_lineages.Rdata")


#set to use all clades
selection="All"
sdat = subAll

# #--- or uncomment chunk below if you want to subset
#pick the clade
selection = "Galaxia"
selection = "Montipora"
selection = "Porites"
selection = "Pocilloporid"

sdat = subAll[!grepl(selection, subAll$b1Clade) & !grepl(selection, subAll$b2Clade),] #for ignoring a clade
sdat = subAll[grepl(selection, subAll$b1Clade) | grepl(selection, subAll$b2Clade),]
head(subAll, n=20)
nrow(dat)
nrow(subAll)

# #-------------

#grab clade options
clades = unique(c(sdat$b1Clade, sdat$b2Clade))
clades
targetClades = clades[!grepl("Anti", clades)]
antiClades = clades[grepl("Anti", clades)]
targetClades
antiClades

#subset dataframes with different associations
tdat = sdat[sdat$b1Clade %in% targetClades & sdat$b2Clade %in% targetClades,]
adat = sdat[sdat$b1Clade %in% antiClades & sdat$b2Clade %in% antiClades,]


#set booleans for which the two clades are different
b1 = sdat$b1Clade %in% targetClades & sdat$b2Clade %in% antiClades
b2 = sdat$b1Clade %in% antiClades & sdat$b2Clade %in% targetClades
b = b1 | b2


#set up bdat
bdat0 = sdat[b,]
bdat0$mod1 = sub("Anti", "", bdat0$b1Clade)
bdat0$mod2 = sub("Anti", "", bdat0$b2Clade)
bdat = bdat0[bdat0$mod1 != bdat0$mod2,]


#check the size of the datasets
#note bdat is larger because of more potential pairs
nrow(bdat)
nrow(tdat)
nrow(adat)
datList = list(tdat, bdat, adat)

#check the pair sets
bpairs = unique(paste(bdat$b1Clade, bdat$b2Clade, sep="_"))
tpairs = unique(paste(tdat$b1Clade, tdat$b2Clade, sep="_"))
apairs = unique(paste(adat$b1Clade, adat$b2Clade, sep="_"))
length(bpairs)
length(tpairs)
length(apairs)

#look at general convergence vs nonconvergence
convergent = unlist(lapply(datList, function(x) sum(x[,'para'])))
not.convergent = unlist(lapply(datList, function(x) sum(x[,'divergent'])))
tbl = rbind(convergent, not.convergent)
labs = c('tt', 'ta', 'aa')
colnames(tbl) = labs

#compare all three
par(mfrow = c(1,3))
barplot(tbl, beside=T, names=labs, legend=F, main=selection)
ratios = (convergent / not.convergent)*100; barplot(ratios, names=labs, ylab='% (Convergent / Non.Convergent)', main = selection)
chisq.test(tbl)

#stats for target-target convergence compared to non-non
tbl2 = tbl[,c('tt', 'aa')]
ratios = (tbl2[1,] / tbl2[2,])*100; barplot(ratios, names=colnames(tbl2), ylab='% (Convergent / Non.Convergent)', main =selection)
fisher.test(tbl2)

#stats for target-target convergence compared to target-non
tbl2 = tbl[,c('tt', 'ta')]
ratios = (tbl2[1,] / tbl2[2,])*100; barplot(ratios, names=colnames(tbl2), ylab='% (Convergent / Non.Convergent)', main =selection)
fisher.test(tbl2)

