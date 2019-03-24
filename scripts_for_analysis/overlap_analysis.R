#overlap_analysis.R
#This is a long version of branch_sites_overlap.R, written first,
# it has extra bits that aren't necessary, but are kept for reference.

#---- LIBS ----#
library(tidyverse)
library(ggplot2)
library(cowplot)
setwd("~/gitreps/convergent_evo_coral")
source('scripts_for_analysis/coral_reproduction_functions.R')



#--- SELECT CUTOFF SET ---#

source("scripts_for_analysis/cutoff_sets/relaxed_cutoff_set.R")          #RELAXED CUTOFFS
source("scripts_for_analysis/cutoff_sets/relaxed_no_brood_cutoff_set.R")   #RELAXED NO BROOD CUTOFFS
source("scripts_for_analysis/cutoff_sets/moderate_cutoff_set.R")           #MODERATE CUTOFFS
source("scripts_for_analysis/cutoff_sets/moderate_no_brood_cutoff_set.R")  #MODERATE NO BROODER PAIRS CUTOFFS
source("scripts_for_analysis/cutoff_sets/semi_strict_cutoff_set.R")        #SEMI-STRICT CUTOFF 
source("scripts_for_analysis/cutoff_sets/strict_cutoff_set.R")             #STRICT CUTOFF




###################################
############ LOAD DATA ############
###################################

#load the three types of data to be compared:
#branch-sites statistics for whole genes (bs)
#flagged sites identified for positive selection (flag)
#convergence data identified from ancestral reconstructions (conv)

#--- LOAD BRANCH-SITES DATA ---#

#load bs data for target clades
gdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_Gacrhelia.tsv", "Galaxia")
mdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_Montipora.tsv", "Montipora")
pdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_Porites.tsv", "Porites")
qdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_Pocillopora.tsv", "Pocilloporid")
vdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_allVertical.tsv", "AllTarget")

#load bs data for anti-target clades
a_gdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_anti_Gacrhelia.tsv", "AntiGalaxia")
a_mdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_anti_Montipora.tsv", "AntiMontipora")
a_pdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_anti_Porites.tsv", "AntiPorites")
a_qdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_anti_Pocillopora.tsv", "AntiPocilloporid")
a_vdat = readin_branch_sites_lrt_results("branch_sites_tests/bs_lrt_anti_allVertical.tsv", "AntiAllTarget")

#build single df
bs_dat = rbind(gdat, mdat, pdat, qdat, vdat, a_gdat, a_mdat, a_pdat, a_qdat, a_vdat) %>% as_tibble
remove(gdat, mdat, pdat, qdat, vdat, a_gdat, a_mdat, a_pdat, a_qdat, a_vdat)
ordered.clades = rev(c('AllTarget', 'AntiAllTarget', 'Montipora', 'AntiMontipora', 'Galaxia', 'AntiGalaxia', 'Porites', 'AntiPorites', 'Pocilloporid', 'AntiPocilloporid'))
bs_dat$clade = factor(bs_dat$clade, levels = ordered.clades)
bs_dat <- bs_dat %>% rename(ortholog = 'contig')


#--- LOAD SITES FLAGGED SITES ---#

#these are the sites flagged for positive selection
#with bayes empirical bayes in the codeml outputs
#build these files with extract_sig_bs_positions.py

#load beb targets
gflag = readin_branch_sites_lrt_results("branch_sites_tests/gachrelia_sigBsSites.tsv", "Galaxia")
mflag = readin_branch_sites_lrt_results("branch_sites_tests/montipora_sigBsSites.tsv", "Montipora")
pflag = readin_branch_sites_lrt_results("branch_sites_tests/porites_sigBsSites.tsv", "Porites")
qflag = readin_branch_sites_lrt_results("branch_sites_tests/pocillopora_sigBsSites.tsv", "Pocilloporid")
vflag = readin_branch_sites_lrt_results("branch_sites_tests/allVertical_sigBsSites.tsv", "AllTarget")

#load beb anti-targets
a_gflag = readin_branch_sites_lrt_results("branch_sites_tests/anti_gachrelia_sigBsSites.tsv", "AntiGalaxia")
a_mflag = readin_branch_sites_lrt_results("branch_sites_tests/anti_montipora_sigBsSites.tsv", "AntiMontipora")
a_pflag = readin_branch_sites_lrt_results("branch_sites_tests/anti_porites_sigBsSites.tsv", "AntiPorites")
a_qflag = readin_branch_sites_lrt_results("branch_sites_tests/anti_pocillopora_sigBsSites.tsv", "AntiPocilloporid")
a_vflag = readin_branch_sites_lrt_results("branch_sites_tests/anti_allVertical_sigBsSites.tsv", "AntiAllTarget")

#build single df
flag_dat = rbind(gflag, mflag, pflag, qflag, vflag, a_gflag, a_mflag, a_pflag, a_qflag, a_vflag) %>% as_tibble
remove(gflag, mflag, pflag, qflag, vflag, a_gflag, a_mflag, a_pflag, a_qflag, a_vflag)
flag_dat$clade = factor(flag_dat$clade, levels = ordered.clades)
flag_dat$file = sub('_ALT.codeml', '', flag_dat$file)
flag_dat <- flag_dat %>% rename('ortholog' = 'file')
flag_dat$site = paste(flag_dat$ortholog, flag_dat$position, sep="_")
flag_dat <- flag_dat %>% mutate(flag_phenotype = if_else(grepl('Anti', clade), 'Horizontal', 'Vertical'))


#--- LOAD CONVERGENCE RESUTLS ---#

#This is the ancestral results file generated with parse_multi_rstV2.py, 
#formatted with prepareConvergenceInput.R
ll = load(ancestral.file); print(ll)
conv_dat$site = paste(conv_dat$ortholog, conv_dat$b1pos, sep="_")
head(conv_dat)

####################################################
############ APPLY FILTERING PARAMETERS ############
####################################################

#--- FOR BRANCH SITE TESTS ---#
#note the calls are made separately for each phenotype,
#because there is a possible filter for the number of lineages within the phenotype

#get genes passing filters for verticals
bs_target_passing_genes <- bs_dat %>%
	filter(!grepl("Anti", clade)) %>%
	filter_bs()

#get genes passing filters for horizontals
bs_anti_passing_genes <- bs_dat %>%
	filter(grepl("Anti", clade)) %>%
	filter_bs()

#assemble dataframe with filter calls
	#important factor assigned:
	#bs_tpass (TRUE if Phenotype is vertical, bs test is significant, 
	#and the number of significant reuslts for the gene was above bs.min.lineage)
	#same for bs_apass

bs_filt <- bs_dat %>% 
	mutate(sig = bs_dat[, bs.ptype] < bs.alpha) %>%
	mutate(bs_phenotype = if_else(grepl('Anti', clade), 'Horizontal','Vertical')) %>% 
	mutate(bs_tpass = 
		bs_phenotype=='Vertical' &
		sig &
		ortholog %in% bs_target_passing_genes) %>% 
	mutate(bs_apass =
		bs_phenotype=='Horizontal' &
		sig &
		ortholog %in% bs_anti_passing_genes) %>% 
	mutate(bs_pass = bs_tpass | bs_apass) %>%
	select(c('ortholog', 'p.values', 'adj.p', 'bs_phenotype', 'clade', 'bs_pass', 'bs_tpass', 'bs_apass')) 
head(bs_filt)

#how many genes pass these filters for the two phenotpyes?
bs_t_genes = unique(bs_filt$ortholog[bs_filt$bs_tpass])
bs_a_genes = unique(bs_filt$ortholog[bs_filt$bs_apass])
length(bs_t_genes)
length(bs_a_genes)



#plot barplot
bs_filt %>% 
	mutate(`Significant` = 
		(bs_phenotype=='Vertical' & bs_tpass) | 
		(bs_phenotype=='Horizontal' & bs_apass) ) %>%
	ggplot() + 
	geom_bar(aes(x=clade, fill=`Significant`)) +
	labs(title='Tests for positive selection', x='Lineage', y='Gene count') +
	coord_flip() +
	theme(axis.text.x= element_text(angle=10))


#--- FOR FLAGGED SITES ---#
#For a flagged site to pass, it must:
#1. Be found in a gene with a branch-site p-value below threshold set as bs.gene.pvalue.cutoff
#2. Must have posterior greater than or equal to value set by bs.site.post.cutoff

#JOIN UP THE FLAG DATA WITH THE BS_PASSABLE COLUMN
#(T/F for whether the bs test for that clade and gene was below bs.gene.pvalue.cutoff)

#make trim df
bs_ftrim <- bs_filt %>% 
	mutate(bs_passable = p.values < bs.gene.pvalue.cutoff) %>%
	select(clade, ortholog, bs_passable)

#join up and call filter passing flagged sites
flag_filt0 <- flag_dat %>%
	left_join(bs_ftrim, by = c('clade', 'ortholog')) %>%
	as_tibble


#IDENTIFY FLAGGED SITES PASSING FILTERS

#for verticals
flag_target_passing_sites <- flag_filt0 %>%
	filter(flag_phenotype=='Vertical') %>%
	filter_flags()


#for horizontals
flag_anti_passing_sites <- flag_filt0 %>%
	filter(flag_phenotype=='Horizontal') %>%
	filter_flags()

#assemble dataframe with filter calls
flag_filt <- flag_filt0 %>% 
	mutate(flag_phenotype = if_else(grepl('Anti', clade), 'Horizontal','Vertical')) %>% 
	mutate(flag_tpass = 
		flag_phenotype =='Vertical' &
	 	posterior >= bs.site.post.cutoff & 
	 	site %in% flag_target_passing_sites) %>%
	mutate(flag_apass = 
		flag_phenotype =='Horizontal' & 
		posterior >= bs.site.post.cutoff & 
		site %in% flag_anti_passing_sites) %>%
	mutate(flag_pass = flag_tpass | flag_apass) %>%
	select(c('ortholog', 'position', 'site', 'aa', 'posterior', 'flag_phenotype', 'clade', 'flag_pass')) %>%
	as_tibble
head(flag_filt)

#how many sites pass filters for the two phenotypes?
sum(flag_filt$flag_phenotype == 'Vertical' & flag_filt$flag_pass)
sum(flag_filt$flag_phenotype == 'Horizontal' & flag_filt$flag_pass)

#how many genes have at least one site passing filters for the two phenotypes?
length(unique(flag_filt$ortholog[flag_filt$flag_pass & flag_filt$flag_phenotype == 'Vertical']))
length(unique(flag_filt$ortholog[flag_filt$flag_pass & flag_filt$flag_phenotype == 'Horizontal']))

#are posteriors as we expect?
flag_filt %>% filter(flag_pass) %>% summarize(minpost=min(posterior)) >= bs.site.post.cutoff

#plot passing rates
flag_filt %>%
	mutate(`Filter Passing`= factor(flag_pass, levels=c(FALSE, TRUE))) %>%
	ggplot() +	geom_bar(aes(x=clade, fill=`Filter Passing`)) + 
	labs(title=paste("Positive sites"), x='Clade', y='Substitution count') +
	coord_flip() +
	theme(axis.text.x= element_text(angle=10))


#--- FOR CONVERGENCE ---#
#For a convergence event to pass, it must be:
#1. parallel or convergent
#2. have posteriors greater than or equal to anc.post.cutoff
#3. must have occured in greater than or equal to min.convergent.changes lineages (min.convergent.changes = the total number of indepdendent convergences WITHIN THE SAME PHENOTYPE PAIRING (VV, VH, or HH))


#FIRST IDENTIFY FILTER-PASSING SITES FOR EACH PHENOTYPE PAIRING

#set up booleans for target/antiTarget clades
t1 = !grepl('Anti', conv_dat$b1Clade)
t2 = !grepl('Anti', conv_dat$b2Clade)
a1 = !t1
a2 = !t2


#filter for target-target convergence events
conv_tt_filt <- conv_dat %>% 
	filter(t1 &t2) %>%
	filter_convergence()


#filter for target-antiTarget convergence events
conv_ta_filt <- conv_dat %>% 
	filter( (t1 & a2) | (a1 & t2) ) %>%
	filter_convergence()

#filter antiTarget-antiTarget convergence events
conv_aa_filt <- conv_dat %>% 
	filter(a1 & a2) %>%
	filter_convergence()


#ASSEMBLE DATAFRAME WITH FILTER CALLS

#set up convergence pairing factor
conv_dat$conv_pair = 'VV'
conv_dat$conv_pair[(t1 & a2) | (a1 & t2)] <- 'VH'
conv_dat$conv_pair[(a1 & a2) | (a1 & a2)] <- 'HH'
conv_dat$conv_pair = factor(conv_dat$conv_pair, levels=c('VV', 'VH', 'HH'))

#modify change type names
type = as.character(conv_dat$sub.type)
type[type=='convergent']<-'convergent (L>V; S>V)'
type[type=='parallel']<-'parallel (L>V; L>V)'
type[type=='divergent']<-'divergent (L>V; L>S)'
type[type=='none']<-'all different (L>V; S>P)'
ordered.type = c('all different (L>V; S>P)', 'divergent (L>V; L>S)', 'parallel (L>V; L>V)', 'convergent (L>V; S>V)')
ordered.type = c('convergent (L>V; S>V)', 'parallel (L>V; L>V)', 'divergent (L>V; L>S)', 'all different (L>V; S>P)')
type = factor(type, levels= ordered.type)
conv_dat$type = type


#assemble dataframe with filter calls
#important designators are conv_pair (which pair of phenotypes) and
#conv_pass (whether the substitution passed filters)
#also filter away sister comparisons in the VH pairings,
#(eg, remove AntiGalaxia-Galaxia comparisons)

conv_filt <- conv_dat %>%
	rename(pos = b1pos) %>%
	mutate(postPass = b1post >= anc.post.cutoff & b2post >= anc.post.cutoff) %>%
	mutate(conv_ttPass = conv_pair=='VV' & genConvergent & postPass & site %in% conv_tt_filt) %>%
	mutate(conv_taPass = conv_pair=='VH' & genConvergent & postPass & site %in% conv_ta_filt) %>%
	mutate(conv_aaPass = conv_pair=='HH' & genConvergent & postPass & site %in% conv_aa_filt) %>%
	mutate(conv_pass = conv_ttPass | conv_taPass | conv_aaPass) %>%
	mutate(sisterComparison = sub('Anti', '', b1Clade)==sub('Anti', '', b2Clade)) %>%
	filter(!sisterComparison) %>% 
	as_tibble
	

#OPTIONAL REMOVAL OF PARTICULAR PAIRS
#This step is for cases when other derived phenotypes than
#those of interest overlap between clades

#remove clade pairs
if (length(remove_conv_clade_pairs) > 0){
	print("Removing the following clade pairings from convergence dataset:")
	for (rmPair in remove_conv_clade_pairs){
		print('---')
		print(rmPair)
		before=nrow(conv_filt)
		conv_filt <- conv_filt %>%
			filter( !(b1Clade==rmPair[1] & b2Clade==rmPair[2]) ) %>%
			filter( !(b1Clade==rmPair[2] & b2Clade==rmPair[1]) )
		after=nrow(conv_filt)
		print(paste(before-after, 'overlapping substititions removed'))
	}
} else{
	print("Not removing any clade pairs")
}


#remove species
if (length(remove_conv_clade_spp_pairs) > 0){
	print("Removing the following clade-species pairings from convergence dataset:")
	for (rmPair in remove_conv_clade_spp_pairs){
		print('---')
		print(rmPair)
		before=nrow(conv_filt)
		conv_filt <- conv_filt %>%
			filter( !(b1Clade==rmPair[1] & grepl(rmPair[2], b2Name)) ) %>%
			filter( !(b2Clade==rmPair[1] & grepl(rmPair[2], b1Name)) )
		after=nrow(conv_filt)
		print(paste(before-after, 'overlapping substititions removed'))
	}
} else{
	print("Not removing any clade-species pairs")
}



#FINAL COLUMN SELECTION FOR conv_filt
conv_filt <- conv_filt %>%
	select('ortholog', 'pos', 'site', 'b1Clade', 'b2Clade', 'b1post', 'b2post', 'conv_pair', 'type', 'genConvergent', 'conv_pass', 'b1Anc', 'b2Anc', 'b1New', 'b2New') %>%
		as_tibble


#------- PLOT FILTER PASSING FREQUENCIES -------#

#BREAKDOWN OF OVERLAPPING SUBSTITUTIONS BY CLADE-PAIR

#build a single convergence pair factor
order_pair = function(x){
	o=x[order(x)]
	r=paste(o[1], o[2], sep='-')
	return(r)
}
cp = apply(conv_filt[,c('b1Clade', 'b2Clade')], 1, function(x) order_pair(x))
table(cp)

#create a 'general' clade pairing vector, where vertical and horiontal 
#pairs have same identity (remove Anti- prefixes)
gcp <- conv_filt %>%
	mutate(gb1=sub('Anti', '', b1Clade), gb2=sub('Anti', '', b2Clade)) %>%
	select(c(gb1, gb2)) %>%
	apply(1, function(x) order_pair(x))
table(gcp)
conv_filt$gcp=gcp


#assign to dataframe and sort to get order by phenotype pairings
conv_filt$cp = cp
ordered_pairs <- conv_filt %>% 
	arrange(desc(conv_pair), cp) %>% 
	select(cp) %>% unique()
conv_filt$cp = factor(conv_filt$cp, levels=ordered_pairs$cp)


#plot all types
conv_filt %>% 
	rename(`Change type` = type) %>%
	ggplot() +
		geom_bar(aes(x=cp, fill=`Change type`)) +
		labs(x='Lineage Pair', y='Overlapping substititions', title='Breakdown of overlapping substitions between lineages') +
		coord_flip()

#plot convergence event frequency
conv_filt %>% 
	mutate(genConvergent = 
		if_else(genConvergent,
			'convergence',
			'no convergence'
			)) %>%
	mutate(genConvergent = factor(genConvergent, levels=c('convergence', 'no convergence'))) %>%
	rename(`Change type`=genConvergent) %>%
	ggplot() +
		geom_bar(aes(x=cp, fill= `Change type`)) +
		labs(x='Lineage Pair', y='Overlapping substititions', title='Frequency of convergence events') +
		coord_flip()



#frequency of genes with overlapping substitions for each pair
conv_filt %>%
	group_by(cp) %>%
	summarize(Ngenes=length(unique(ortholog))) %>%
	ggplot() +
		geom_bar(aes(x=cp, y=Ngenes), stat='identity') +
		labs(x='Lineage Pair', y='Gene count', title='Frequency of genes with convergence events') +
		coord_flip()

#plot frequency of genes with convergence events
conv_filt %>% 
	group_by(cp, conv_pass) %>%
	summarize(Ngenes=length(unique(ortholog))) %>%
	mutate(`Convergence detected`=conv_pass) %>%
	ggplot() +
		geom_bar(aes(x=cp, y=Ngenes, fill= `Convergence detected`), stat='identity') +
		labs(x='Lineage Pair', y='Gene count', title='Frequency of genes with convergence events') +
		coord_flip()




#plot frequency of filter-passing to all convergent subs
conv_filt %>%
	filter(genConvergent) %>%
	rename(`Filter-passing`=conv_pass) %>%
	ggplot() +
		geom_bar(aes(x=cp, fill=`Filter-passing`)) + 
		labs(x='Convergence pair', y='Convergence events', title='Frequency of filter-passing convergence events') +
		coord_flip()

#plot frequency of filter-passing to all convergent subs by phenotype pairing
conv_filt %>%
	filter(genConvergent) %>%
	rename(`Filter-passing`=conv_pass) %>%
	ggplot() +
		geom_bar(aes(x=conv_pair, fill=`Filter-passing`)) + 
		labs(x='Convergence pair', y='Convergence events', title='Frequency of filter-passing convergence events')


#gather percentages of convergence passing subs
pctdf <-conv_filt %>% 
	mutate(conv_pass = 
		if_else(conv_pass,
			'convergence',
			'no convergence'
			)) %>%
	mutate(conv_pass = factor(conv_pass, levels=c('convergence', 'no convergence'))) %>%
	rename(`Change type`= conv_pass) %>%
	group_by(cp, conv_pair, `Change type`) %>%
	summarize(N=n()) %>% 
	mutate(pct_conv = N/sum(N)*100) %>%
	filter(`Change type`=='convergence') %>%
	# rename(`Phenotype Pair` = conv_pair) %>%
	arrange(conv_pair)

#plot by clade pairing
	pctdf %>%
	rename(`Phenotype Pair` = conv_pair) %>%
	ggplot() +
		geom_bar(aes(x=cp, y= pct_conv, fill=`Phenotype Pair`), stat='identity') +
		labs(x='Lineage Pair', y='Percent convergence passing', title='convergent : overlapping substititions') +
		coord_flip()


#plot mean percentage by phenotype pairing
pctdf %>% 
	rename(`Phenotype Pair` = conv_pair) %>%
	rename(`Percent Convergent` = pct_conv) %>%
	ggplot() +
	geom_boxplot(aes(x=`Phenotype Pair`, y=`Percent Convergent`, fill=`Phenotype Pair`))
	
#basic stats
fit <- aov(pct_conv ~ conv_pair, data=pctdf)
summary(fit)


#paired comparison
pair <- conv_filt %>% 
	mutate(conv_pass = 
		if_else(conv_pass,
			'convergence',
			'no convergence'
			)) %>%
	mutate(conv_pass = factor(conv_pass, levels=c('convergence', 'no convergence'))) %>%
	rename(`Change type`= conv_pass) %>%
	group_by(gcp, conv_pair, `Change type`) %>%
	summarize(N=n()) %>% 
	mutate(pct_conv = N/sum(N)*100) %>%
	filter(`Change type`=='convergence') %>%
	# rename(`Phenotype Pair` = conv_pair) %>%
	arrange(conv_pair)

#compare VV and HH
tin <- pair[,c('conv_pair', 'pct_conv', 'gcp')] %>%
	filter(conv_pair!='VH') %>%
	spread(conv_pair, pct_conv) 
t.test(x=tin$VV, y=tin$HH, paired=T)

#compare VV and VH
tin <- pair[,c('conv_pair', 'pct_conv', 'gcp')] %>%
	filter(conv_pair!='HH') %>%
	spread(conv_pair, pct_conv) 
t.test(x=tin$VV, y=tin$VH, paired=T)

#absolute counts divided by sub type
c1<-conv_filt %>% 
	rename(`Change type` = type) %>%
	ggplot() +
	geom_bar(aes(x=conv_pair, fill=`Change type`)) + 
	labs(x='Convergence pairing', y='Changes at same position', title='Overlapping substitutions by type', subtitle='(All overlapping substititions)')
	
#absolute counts divided by sub type for subs passing posterior filter
subtitle = paste(paste('(substitution posterior >=', anc.post.cutoff), ')', sep='')
c2<-conv_filt %>% 
	filter(b1post & b2post >= anc.post.cutoff) %>%
	rename(`Change type` = type) %>%
	ggplot() +
	geom_bar(aes(x=conv_pair, fill=`Change type`)) + 
	labs(x='Convergence pairing', y='Changes at same position', title='Overlapping substitutions by type', subtitle=subtitle)
	
plot_grid(c1, c2)



#absolute counts convergent/not convergence for all subs
c1<-conv_filt %>% 
	mutate(genConvergent = 
		if_else(genConvergent,
			'convergence',
			'no convergence'
			)) %>%
	mutate(genConvergent = factor(genConvergent, levels=c('convergence', 'no convergence'))) %>%
	rename(`Change type`=genConvergent) %>%
	ggplot() +
	geom_bar(aes(x=conv_pair, fill= `Change type`)) + 
	labs(x='Phenotype pairing', y='Changes at same position', title='Frequency of convergence events', subtitle='(All Changes)')

#absolute counts convergent filter passing or not passing
c2<-conv_filt %>% 
	mutate(`Convergence` = 
		if_else(conv_pass,
			'pass',
			'fail'
			)) %>%
	mutate(`Convergence` = factor(`Convergence`, levels=c('pass', 'fail'))) %>%
	ggplot() +
	geom_bar(aes(x=conv_pair, fill= `Convergence`)) + 
	labs(x='Phenotype pairing', y='Changes at same position', title='Frequency of convergence events', subtitle='(Passing convergence parameters)')

plot_grid(c2)

#plot proportions for convergence frequency
conv_filt %>% 
	group_by(conv_pair, genConvergent) %>%
	summarize(N=n()) %>%
	mutate(`Percent Convergence`=N/sum(N)*100) %>%
	filter(genConvergent) %>% 
	ggplot() +
		geom_bar(aes(x=conv_pair, y=`Percent Convergence`, fill=conv_pair), stat='identity') +
		labs(y='Convergence (% overlapping changes)', title='Convergence among overlapping substititions', subtitle='(All substititions)')
		

#plot proportions for convergence filter passing
conv_filt %>% 
	group_by(conv_pair, conv_pass) %>%
	summarize(N=n()) %>%
	mutate(`Percent Convergence`=N/sum(N)*100) %>%
	filter(conv_pass) %>% 
	ggplot() +
		geom_bar(aes(x=conv_pair, y=`Percent Convergence`, fill=conv_pair), stat='identity') +
		labs(y='Convergence (% overlapping changes)', title='Convergence among overlapping substititions', subtitle='(Filter passing substititions)', x='Convergence pair')


#stats on proportion
ft_in <- conv_filt %>%
	filter(conv_pair != 'VH') %>%
	mutate(conv_pair = factor(conv_pair, levels=c('VV', 'HH'))) %>%
	mutate(convergent = factor(conv_pass, levels=c(TRUE, FALSE))) %>%
	select(c('conv_pair', 'convergent'))
fisher.test(ft_in$conv_pair, ft_in$convergent)

#prove to myselft this works
x<- conv_filt %>% 
	filter(conv_pair != 'VH') %>%
	mutate(conv_pair = factor(conv_pair, levels=c('VV', 'HH'))) %>%
	mutate(convergent = factor(conv_pass, levels=c(TRUE, FALSE))) %>%
	group_by(conv_pair, convergent) %>%
	summarize(N=n())

t=matrix(x$N, nrow=2)
fisher.test(t)

x=table(paste(conv_filt$conv_pair, conv_filt$conv_pass, sep='_'))
vv=x[c('VV_TRUE', 'VV_FALSE')]
hh=x[c('HH_TRUE', 'HH_FALSE')]
t=as.table(rbind(vv,hh))
rownames(t) = c('VV', 'HH')
colnames(t) = c('Pass', 'Fail')
t
fisher.test(t)



################################################
##### BRANCH SITES AND CONVERGENCE OVERLAP #####
################################################
#add columns to the conv_filt dataframe for clade1 and clade2
#of the convergence pair that indicates if
#the site if found in a gene that passed filters in 
#the bs data for the clade itself, or the appropriate 'all' clade

#relevent dataframes:
head(bs_filt)
head(conv_filt)


#Build a list of bs filter passing genes for each clade
clades = unique(bs_filt$clade)
posGeneList = lapply(clades, function(x) gather_pos_genes(x))
names(posGeneList) = clades

#double-check things make sense
length(posGeneList[['Galaxia']])
bs_filt %>% 
	filter(clade=='Galaxia') %>%
	filter(bs_pass) %>%
	nrow()
length(unique(bs_filt$ortholog[bs_filt$clade=='Galaxia' & bs_filt$bs_pass]))

#use cross_check_clade1() and cross_check_clade2() to assign if 
#subs are also found in positively selected genes for the given clade or 
#the appropriate 'all' group for the phenotype
#spot checked this on 10-13-18
conv_filt$c1PosGene = apply(conv_filt, 1, function(x) cross_check_bs_clade1(x))
conv_filt$c2PosGene = apply(conv_filt, 1, function(x) cross_check_bs_clade2(x))


#PLOT FREQUENCY FOR VERTICALS ONLY

#counts of convergence events (substitions)
ordered <- conv_filt %>%
	group_by(cp, conv_pass) %>%
	summarize(N=n()) %>%
	filter(conv_pass) %>%
	arrange(N)
	
#plot frequency of substitition and bs passing
conv_filt %>%
	filter(conv_pair=='VV') %>%
	filter(conv_pass) %>%
	mutate(posin = as.character(c1PosGene + c2PosGene)) %>%
	mutate(posin = if_else(posin=='0', 'none', posin)) %>% 
	mutate(posin = if_else(posin=='1', 'one', posin)) %>% 
	mutate(posin = if_else(posin=='2', 'both', posin)) %>% 
	mutate(posin = factor(posin, levels=c('none', 'one', 'both'))) %>%
	rename(`Positive selection in`=posin) %>%
	mutate(cp = factor(cp, levels=ordered$cp)) %>%
	ggplot() +
		geom_bar(aes(x=cp, fill=`Positive selection in`)) +
		labs(x='Convergence pair', y='Substitition count', title='Convergence and positive selection') +
		coord_flip()

#counts for genes with at least one convergence event
ordered <- conv_filt %>%
	filter(conv_pair=='VV') %>%
	group_by(cp, conv_pass) %>%
	summarize(Ngenes=length(unique(ortholog))) %>%
	filter(conv_pass) %>%
	arrange(Ngenes)

	
conv_filt %>%
	filter(conv_pair=='VV') %>%
	filter(conv_pass) %>%
	mutate(posin = as.character(c1PosGene + c2PosGene)) %>%
	mutate(posin = if_else(posin=='0', 'none', posin)) %>% 
	mutate(posin = if_else(posin=='1', 'one', posin)) %>% 
	mutate(posin = if_else(posin=='2', 'both', posin)) %>% 
	mutate(posin = factor(posin, levels=c('none', 'one', 'both'))) %>%
	mutate(cp = factor(cp, levels=ordered$cp)) %>%
	group_by(cp, posin) %>% 
	summarize(Ngenes=length(unique(ortholog))) %>% 
	rename(`Positive selection in`=posin) %>% 
	ggplot() +
		geom_bar(aes(x=cp, y=Ngenes, fill=`Positive selection in`), stat='identity') +
		labs(x='Convergence pair', y='Genes with convergence events', title='Convergence and positive genes') +
		coord_flip()

#COMPARE CO-OCCURANCE AMONG ALL TYPES

#plot frequencies for all (substitutions)
conv_filt %>%
	filter(conv_pass) %>%
	mutate(`Positive tests` = factor(c1PosGene + c2PosGene)) %>%
	ggplot() +
		geom_bar(aes(x=cp, fill=`Positive tests`)) +
		labs(x='Convergence pair', y='Convergence events', title='Convergence and positive selection') +
		coord_flip()


#plot frequencies for all (gene counts)
conv_filt %>%
	filter(conv_pass) %>%
	mutate(`Positive tests` = factor(c1PosGene + c2PosGene)) %>%
	group_by(cp, `Positive tests`) %>%
	summarize(Ngenes=length(unique(ortholog))) %>%
	ggplot() +
		geom_bar(aes(x=cp, y= Ngenes, fill=`Positive tests`), stat='identity') +
		labs(x='Convergence pair', y='Genes with convergence', title='Convergence and positive selection') +
		coord_flip()

#plot percentages of positive selection among convergent genes
conv_filt %>%
	filter(conv_pass) %>%
	mutate(`Positive tests` = factor( (c1PosGene + c2PosGene)>0 )) %>%
	group_by(cp, conv_pair, `Positive tests`) %>%
	summarize(Ngene=length(unique(ortholog))) %>% mutate(pctBoth = Ngene / sum(Ngene)*100) %>% 
	filter(`Positive tests`==TRUE) %>%
	rename(`Phenotype pair`= conv_pair) %>%
	ggplot() +
		geom_bar(aes(x=cp, y=pctBoth, fill=`Phenotype pair`), stat='identity') +
		labs(x='Convergence pair', y='Percentage of convergent genes', title='Positive selection among convergent genes') +
		coord_flip()



#set up percentages for boxplot (genes)
pct<-conv_filt %>%
	filter(conv_pass) %>%
	mutate(`Positive tests` = factor( (c1PosGene + c2PosGene)>0 )) %>%
	group_by(cp, conv_pair, `Positive tests`) %>%
	summarize(Ngene=length(unique(ortholog))) %>% mutate(pctBoth = Ngene / sum(Ngene)*100) %>%
	filter(`Positive tests`==TRUE) %>%
	rename(`Phenotype pair`= conv_pair)

#plot boxplot
pct %>%
	ggplot() +
		geom_boxplot(aes(x=`Phenotype pair`, y=pctBoth, fill=`Phenotype pair`)) +
		labs(x='Phenotype pair', y='Percentage of convergent genes', title='Positive selection among convergent genes')

#stats on this
fit<-aov(pctBoth~`Phenotype pair`, data=pct)
summary(fit)


#PLOT PERCENTAGE CONVERGENCE AND POSITIVE AMONG ALL GENES TESTED FOR EACH PAIR

#first get the counts of genes that show both
Nboth <- conv_filt %>%
	filter(conv_pass) %>%
	mutate(`Positive tests` = factor( (c1PosGene + c2PosGene)>0 )) %>%
	group_by(cp, conv_pair, `Positive tests`) %>% 
	summarize(Ngene=length(unique(ortholog))) %>%
	filter(`Positive tests`==TRUE)
Nboth %>% data.frame()

#now get the number of genes with branch-site results from both
#(this is the number of genes that could potentially have convergence,
#since both the lineages were present)
pair_arr=c()
ngene_arr=c()
for (c1 in clades){
	for (c2 in clades){
		pair=paste(c1, c2, sep='-')
		sub1 <- bs_filt %>%
			filter(clade == c1)
		sub2 <- bs_filt %>%
			filter(clade == c2)
		sub <- sub1[sub1$ortholog %in% sub2$ortholog,]
		ngenes = length(unique(sub$ortholog))
		pair_arr=append(pair_arr, pair)
		ngene_arr=append(ngene_arr, ngenes)
	}
}
Ngenes = data.frame('cp'=pair_arr, 'nPossible'=ngene_arr)
mdat <- merge(Ngenes, Nboth, by = 'cp') %>%
	mutate(pct=Ngene/nPossible*100) %>%
	arrange(desc(conv_pair),pct) %>%
	mutate(cp=factor(cp, levels=cp))

#plot percentage by convergence pairing
mdat %>%
	mutate(pct = Ngene / nPossible*100) %>%
	rename(`Phenotype pair`=conv_pair) %>%
	ggplot() +
		geom_bar(aes(x=cp, y=pct, fill=`Phenotype pair`), stat='identity') +
		labs(x='Convergence pair', y='Perecentage of tested genes', title='Frequency of convergence and positive selection') +
		coord_flip()

#plot percentage boxplots
mdat %>%
	mutate(pct = Ngene / nPossible*100) %>%
	rename(`Phenotype pair`=conv_pair) %>%
	ggplot() +
		geom_boxplot(aes(x=`Phenotype pair`, y=pct, fill=`Phenotype pair`)) +
		labs(x='Phenotype pair', y='Perecentage of tested genes', title='Frequency of convergence and positive selection')

#stats on this
fit <- aov(pct ~ conv_pair, data=mdat)
summary(fit)

#t-test for VV vs HH
vv=mdat %>% filter(conv_pair=='VV')
hh=mdat %>% filter(conv_pair=='HH')
t.test(x=vv$pct, y=hh$pct)

##################################################
##### FLAGGEDE SITES AND CONVERGENCE OVERLAP #####
##################################################

#relevent dataframes:
head(flag_filt)
head(conv_filt)

#Build a list of flag filter passing sites for each clade
clades = unique(flag_filt$clade)
gather_flag_sites = function(x){
	csub <- flag_filt %>%
		filter(clade==x & flag_pass)
	return(unique(csub$site))
}
flagSiteList = lapply(clades, function(x) gather_flag_sites(x))
names(flagSiteList) = clades


#double-check things make sense
cld = 'AllTarget'
head(flagSiteList[[cld]])
length(flagSiteList[[cld]])
flag_filt %>% 
	filter(clade==cld) %>%
	filter(flag_pass) %>%
	nrow()
length(unique(flag_filt$site[flag_filt$clade==cld & flag_filt$flag_pass]))

# #use cross_check_clade1() and cross_check_clade2() to assign if 
# # subs are also flagged for positive selection in that clade
# # this takes a while
conv_filt$c1FlagSite = apply(conv_filt, 1, function(x) cross_check_flags_clade1(x))
conv_filt$c2FlagSite = apply(conv_filt, 1, function(x) cross_check_flags_clade2(x))


#load the flag calls
fileName = paste(results.path, 'conv_filt_flagged.Rdata', sep='/')
save(conv_filt, file=fileName)
load(fileName)

#PLOT FREQUENCY FOR VERTICALS ONLY

#counts of convergence events AND positive site (substitions)
ordered <- conv_filt %>%
	group_by(cp, conv_pass) %>%
	summarize(N=n()) %>%
	filter(conv_pass) %>%
	arrange(N)
	
conv_filt %>%
	filter(conv_pair=='VV') %>%
	filter(conv_pass) %>%
	mutate(posin = as.character(c1FlagSite + c2FlagSite)) %>%
	mutate(posin = if_else(posin=='0', 'none', posin)) %>% 
	mutate(posin = if_else(posin=='1', 'one', posin)) %>% 
	mutate(posin = if_else(posin=='2', 'both', posin)) %>% 
	mutate(posin = factor(posin, levels=c('none', 'one', 'both'))) %>%
	rename(`Positive site in`=posin) %>%
	mutate(cp = factor(cp, levels=ordered$cp)) %>%
	ggplot() +
		geom_bar(aes(x=cp, fill=`Positive site in`)) +
		labs(x='Convergence pair', y='Substitition count', title='Convergence and positive site') +
		coord_flip()


#---- COMPARE CO-OCCURANCE OF CONVERGENCE AND FLAGGING AMONG ALL TYPES BY SUBSTITION ----#

#plot frequencies for convergence events and overlap with flagged sites (substitutions)
conv_filt %>%
	filter(conv_pass) %>%
	mutate(`Positive site` = (c1FlagSite + c2FlagSite)>0 ) %>%
	ggplot() +
		geom_bar(aes(x=cp, fill=`Positive site`)) +
		labs(x='Convergence pair', y='Convergence events', title='Convergence and positive site (substitions)') +
		coord_flip()

#PLOT PERCENTAGE OF CONVERGENT AND FLAGGED SITES (substitutions)
#organize percentages
pct <- conv_filt %>%
	filter(conv_pass) %>%
	mutate(`Positive site` = (c1FlagSite + c2FlagSite)>0 ) %>%
	group_by(cp, conv_pair, `Positive site`) %>%
	summarize(N=n()) %>% mutate(Percentage=N/sum(N)*100) %>%
	filter(`Positive site`==TRUE) %>%
	arrange(desc(conv_pair), Percentage)
pct$cp = factor(pct$cp, levels=pct$cp)


#plot bars	
pct %>%
	rename(`Phenotype pair` = conv_pair) %>%
	ggplot() +
		geom_bar(aes(x=cp, y=Percentage, fill=`Phenotype pair`), stat='identity') +
		labs(x='Convergence pair', y='Percentage of convergence events', title='Convergence and positive site (substitions)') +
		coord_flip()

#plot boxplots
pct %>%
	rename(`Phenotype pair` = conv_pair) %>%
	ggplot() +
		geom_boxplot(aes(x=`Phenotype pair`, y=Percentage, fill=`Phenotype pair`)) +
		labs(x='Convergence pair', y='Percentage of convergence events also positive', title='Convergence and positive site (substitions)')

#---- COMPARE CO-OCCURANCE OF CONVERGENCE AND FLAGGING AMONG ALL TYPES BY GENE ----#

pairs = unique(conv_filt$cp)
gather_convflag_genes = function(x){
	csub <- conv_filt %>%
		mutate(posPass=c1FlagSite | c2FlagSite) %>%
		filter(cp==x & posPass & conv_pass)
	return(unique(csub$ortholog))
}
flagConvGeneList = lapply(pairs, function(x) gather_convflag_genes(x))
names(flagConvGeneList) = pairs


#double-check things make sense
p = 'Galaxia-Montipora'
head(flagConvGeneList[[p]])
length(flagConvGeneList[[p]])
conv_filt %>% 
	filter(cp==p) %>%
	filter(conv_pass & (c1FlagSite | c2FlagSite)) %>%
	summarize(length(unique(ortholog)))


# #use cross_check_convflags() to assign if a gene is among those with both results for a pair
conv_filt$convFlagGene = apply(conv_filt, 1, function(x) cross_check_convflags(x))



#plot frequencies for both results among all genes 
conv_filt %>% 
	group_by(cp, convFlagGene) %>%
	summarize(count=length(unique(ortholog))) %>% 
	mutate(`Conv. + Pos. Gene`=convFlagGene) %>%
	ggplot() +
		geom_bar(aes(x=cp, y=count, fill=`Conv. + Pos. Gene`), stat='identity') +
		coord_flip()



#plot frequencies for both results among convergence genes 
conv_filt %>%
	filter(conv_pass) %>%
	mutate(`Positive site` = (c1FlagSite + c2FlagSite)>0  ) %>%
	group_by(cp, `Positive site`) %>%
	summarize(Ngenes=length(unique(ortholog))) %>%
	ggplot() +
		geom_bar(aes(x=cp, y= Ngenes, fill=`Positive site`), stat='identity') +
		labs(x='Convergence pair', y='Genes with convergence and positive sites', title='Convergence and positive site (genes)') +
		coord_flip()

#plot percentages of convergent genes and convergent genes with flagged sites
pct<-conv_filt %>%
	filter(conv_pass) %>%
	mutate(`Positive site` = (c1FlagSite + c2FlagSite)>0 ) %>%
	group_by(cp, conv_pair, `Positive site`) %>%
	summarize(Ngene=length(unique(ortholog))) %>% 
  mutate(pctBoth = Ngene / sum(Ngene)*100) %>%
	filter(`Positive site`==TRUE) %>%
	rename(`Phenotype pair`= conv_pair) %>%
	arrange(desc(`Phenotype pair`), pctBoth)
pct$cp = factor(pct$cp, levels=pct$cp)


#plot percentage of convergent genes for each pair
pct %>%
	ggplot() +
		geom_bar(aes(x=cp, y=pctBoth, fill=`Phenotype pair`), stat='identity') +
		labs(x='Convergence pair', y='Percentage of convergent genes', title='Genes with convergence and positive sites') +
		coord_flip()

#plot boxplot (genes)
pct %>%
	ggplot() +
		geom_boxplot(aes(x=`Phenotype pair`, y=pctBoth, fill=`Phenotype pair`)) +
		labs(x='Phenotype pair', y='Percentage of convergent genes', title='Convergence and positive site (genes)')

#stats on this
fit<-aov(pctBoth~`Phenotype pair`, data=pct)
summary(fit)


#######################################################
##### FLAGGEDE SITES AND CONVERGENCE AND POSITIVE #####
#######################################################

#set up factor for overlaps of different data types
overlap <- conv_filt %>%
	mutate(positive = ( (c1PosGene & c1FlagSite) | (c2PosGene & c2FlagSite) )) %>% 
	mutate(fullEv = positive & conv_pass) %>% 
	mutate(allEv = 'not convergent or positive') %>% 
	mutate(allEv = if_else(!positive & conv_pass, 'convergent', allEv)) %>%
	mutate(allEv = if_else(positive & !conv_pass, 'positive', allEv)) %>%
	mutate(allEv = if_else(fullEv, 'convergent and positive', allEv)) %>%
	mutate(allEv = factor(allEv, levels=c('not convergent or positive', 'convergent', 'positive', 'convergent and positive')))


#PLOT FREQUENCY FOR VERTICALS ONLY

#counts of convergence events AND positive site (substitions)
ordered <- conv_filt %>%
	filter(conv_pair=='VV') %>%
	group_by(cp) %>%
	summarize(N=n()) %>%
	arrange(N)

#plot absolute for all overlapping substitions
overlap %>%
	filter(conv_pair=='VV') %>%
	rename(`Overlapping substitions:`= allEv) %>% 
	mutate(cp = factor(cp, levels=ordered$cp)) %>%
	ggplot() +
		geom_bar(aes(x=cp, fill=`Overlapping substitions:`)) +
		labs(x='Convergence pair', y='Overlapping substitions', title='Convergence & positive selection among sites') +
		coord_flip()


#---- LOOK AT BREAKDOWN BY SUBSTITION ----#

#plot breakdown for all pairs for all overlapping subs
overlap %>%
	rename(`Overlapping substitions:`= allEv) %>% 
	ggplot() +
		geom_bar(aes(x=cp, fill=`Overlapping substitions:`)) +
		labs(x='Convergence pair', y='Overlapping substitions', title='Convergence & positive selection among overlapping substitions') +
		coord_flip()


#plot breakdown for all pairs for all overlapping subs without 'nones'
overlap %>%
	filter(allEv!='not convergent or positive') %>%
	rename(`Overlapping substitions:`= allEv) %>% 
	ggplot() +
		geom_bar(aes(x=cp, fill=`Overlapping substitions:`)) +
		labs(x='Convergence pair', y='Overlapping substitions', title='Convergence & positive selection among overlapping substitions') +
		coord_flip()


#plot percentage both convergent and positive
over_pct <- overlap %>%
	group_by(cp, conv_pair, allEv) %>%
	summarize(N=n()) %>% mutate(Percentage=N/sum(N)*100) %>%
	filter(allEv=='convergent and positive') %>%
	arrange(desc(conv_pair), Percentage) %>%
	data.frame() %>%
	mutate(cp=factor(cp, levels=cp)) %>%
	rename(`Phenotype pair` = conv_pair)

#plot percentage by pair
over_pct %>%
	ggplot() +
		geom_bar(aes(x=cp, y=Percentage, fill=`Phenotype pair`), stat='identity') +
		labs(x='Convergence pair', y='Overlapping substitions', title='Convergence & positive selection among overlapping substitions') +
		coord_flip()


#boxplot by phenotype pair
over_pct %>%
	ggplot() +
	geom_boxplot(aes(x=`Phenotype pair`, y=Percentage, fill=`Phenotype pair`))



#---- LOOK AT BREAKDOWN BY GENE ----#

#make lists of genes for each pairing that have all three pieces of evidence
gather_allEv_genes = function(x){
	csub <- overlap %>%
		filter(cp==x & allEv=="convergent and positive")
	return(unique(csub$ortholog))
}

pairs = unique(overlap$cp)
allEvList = lapply(pairs, function(x) gather_allEv_genes(x))
names(allEvList) = pairs
head(allEvList[['Galaxia-Porites']])
length(allEvList[['Galaxia-Porites']])


#assign if genes are flagged for all evidence types in any substition
overlap$geneAll = apply(overlap, 1, function(x) cross_check_all(x))

#double-check
length(unique(overlap$ortholog[overlap$geneAll & overlap$cp=='Galaxia-Porites']))


#plot breakdown for all pairs for all overlapping subs by gene
overlap %>%
	group_by(cp, conv_pair, geneAll) %>%
	summarize(N=length(unique(ortholog))) %>%
	rename(`Convergent and positive`=geneAll) %>%
	ggplot() +
		geom_bar(aes(x=cp, y=N, fill=`Convergent and positive`), stat='identity') +
		labs(x='Convergence pair', y='Overlapping substitions', title='Convergence & positive selection among overlapping substitions') +
		coord_flip()


#plot percentage both convergent and positive
overall_pct <- overlap %>%
	group_by(cp, conv_pair, geneAll) %>%
	summarize(N=length(unique(ortholog))) %>%
	mutate(Percentage=N/sum(N)*100) %>%
	filter(geneAll) %>%
	arrange(desc(conv_pair), Percentage) %>%
	data.frame() %>%
	mutate(cp=factor(cp, levels=cp)) %>%
	rename(`Phenotype pair` = conv_pair)


#plot percentages by convergence pair
overall_pct %>%
	ggplot() +
		geom_bar(aes(x=cp, y=Percentage, fill=`Phenotype pair`), stat='identity') +
		labs(x='Convergence pair', y='Percentage of tested genes', title='Genes with convergence & positive selection') +
		coord_flip()
		
		
#plot boxplot
overall_pct %>%
	ggplot() +
	geom_boxplot(aes(x=`Phenotype pair`, y=Percentage, fill=`Phenotype pair`))

#run t-test comparing VV and HH
vv=overall_pct[overall_pct[,2]=='VV', 'Percentage']
hh=overall_pct[overall_pct[,2]=='HH', 'Percentage']
t.test(vv, hh, alternative='greater')



############################################
######### OUTPUT FOR GO ENRICHMENT #########
############################################
#how many genes had tests for each phenotype?
#isolate the total number of genes with bs data for each clade type
target_bs_genes <- bs_dat %>% 
	mutate(target = !grepl('Anti', clade)) %>%
	filter(target) %>% 
	select(ortholog) %>%
	unique()
anti_bs_genes <- bs_dat %>% 
	mutate(target = grepl('Anti', clade)) %>%
	filter(target) %>%
	select(ortholog) %>%
	unique()
all_target_genes = target_bs_genes$ortholog
all_anti_genes = anti_bs_genes$ortholog
all_ta_genes = unique(append(all_target_genes, all_anti_genes))

#grab the counts of genes tested for each
nt = length(all_target_genes)
na = length(all_anti_genes)
nt
na


setwd("~/gitreps/convergent_evo_coral")
head(conv_filt)
length(all_ta_genes)
length(all_anti_genes)
length(all_target_genes)

#Select the set you want
setwd("~/gitreps/convergent_evo_coral")
select_pair_set = c('VV', 'HH')
select_phenotype_set = c('Vertical', 'Horizontal')
allGeneListPair = list('VV' = all_target_genes, 'HH' = all_anti_genes)
allGeneListPhenotype = list('Vertical' = all_target_genes, 'Horizontal' = all_anti_genes)

gids = read.table("ortholog_tables/singleCopyAnnotationsWithDescriptions.tsv", sep='\t', quote="", header = T)

i=1
for (i in seq(1, length(select_pair_set))){
	select_pair = select_pair_set[i]
	select_phenotype = select_phenotype_set[i]


	#positive
	pc <-bs_filt %>% 
		filter(bs_phenotype == select_phenotype) %>%
		filter(bs_pass) %>%
		select(ortholog)
	pc_genes = unique(pc$ortholog)
	allGenes = allGeneListPhenotype[[select_phenotype]]
	not_pc_genes = allGenes[! allGenes %in% pc_genes]
	pc_score = rep(1, length(pc_genes))
	not_score = rep(0, length(not_pc_genes))
	goin = data.frame('gene' = c(pc_genes, not_pc_genes), 'score' = c(pc_score, not_score))
	head(goin)
	table(goin$score)
	outName = paste(results.path, paste(paste('GO_MWU/GOMWU_input1_branchSites', select_phenotype, sep=''), 'csv', sep='.'), sep='/')
	write.csv(goin, file=outName, quote=F, row.names=F)
	gene_names = gids %>% filter(ortho %in% pc_genes)
	goutName = paste(results.path, paste(select_pair, 'branch_site_passing_genes.tsv', sep='_'), sep=paste('/'))
	write.table(gene_names, file=goutName, sep='\t', quote=F)
	
	#covnergent
	pc <-conv_filt %>% 
		filter(conv_pair== select_pair) %>%
		filter(conv_pass) %>%
		select(ortholog)
	pc_genes = unique(pc$ortholog)
	allGenes = allGeneListPair[[select_pair]]
	not_pc_genes = allGenes[! allGenes %in% pc_genes]
	pc_score = rep(1, length(pc_genes))
	not_score = rep(0, length(not_pc_genes))
	goin = data.frame('gene' = c(pc_genes, not_pc_genes), 'score' = c(pc_score, not_score))
	head(goin)
	table(goin$score)
	outName = paste(results.path, paste(paste('GO_MWU/GOMWU_input2_covergence', select_pair, sep=''), 'csv', sep='.'), sep='/')
	write.csv(goin, file=outName, quote=F, row.names=F)
	gene_names = gids %>% filter(ortho %in% pc_genes)
	goutName = paste(results.path, paste(select_pair, 'convergence_passing_genes.tsv', sep='_'), sep=paste('/'))
	write.table(gene_names, file=goutName, sep='\t', quote=F, row.names=F)
	
	
	#positive and convergent
	pc <-conv_filt %>% 
		filter(conv_pair== select_pair) %>%
		filter(conv_pass & (c1PosGene | c2PosGene) ) %>%
		select(ortholog)
	pc_genes = unique(pc$ortholog)
	allGenes = allGeneListPair[[select_pair]]
	not_pc_genes = allGenes[! allGenes %in% pc_genes]
	pc_score = rep(1, length(pc_genes))
	not_score = rep(0, length(not_pc_genes))
	goin = data.frame('gene' = c(pc_genes, not_pc_genes), 'score' = c(pc_score, not_score))
	head(goin)
	table(goin$score)
	outName = paste(results.path, paste(paste('GO_MWU/GOMWU_input3_bs_convergence_overlap', select_pair, sep=''), 'csv', sep='.'), sep='/')
	write.csv(goin, file=outName, quote=F, row.names=F)
	gene_names = gids %>% filter(ortho %in% pc_genes)
	goutName = paste(results.path, paste(select_pair, 'branch_site_and_convergence_passing_genes.tsv', sep='_'), sep=paste('/'))
	write.table(gene_names, file=goutName, sep='\t', quote=F, row.names=F)
	
	#flagged and convergent
	pc <- conv_filt %>% 
		filter(conv_pair== select_pair) %>%
		filter(convFlagGene) %>%
		select(ortholog)
	pc_genes = unique(pc$ortholog)
	allGenes = allGeneListPair[[select_pair]]
	not_pc_genes = allGenes[! allGenes %in% pc_genes]
	pc_score = rep(1, length(pc_genes))
	not_score = rep(0, length(not_pc_genes))
	goin = data.frame('gene' = c(pc_genes, not_pc_genes), 'score' = c(pc_score, not_score))
	head(goin)
	table(goin$score)
	outName = paste(results.path, paste(paste('GO_MWU/GOMWU_input4_flagged_convergence_overlap', select_pair, sep=''), 'csv', sep='.'), sep='/')
	write.csv(goin, file=outName, quote=F, row.names=F)
	gene_names = gids %>% filter(ortho %in% pc_genes)
	goutName = paste(results.path, paste(select_pair, 'flagged_and_convergence_passing_genes.tsv', sep='_'), sep=paste('/'))
	write.table(gene_names, file=goutName, sep='\t', quote=F, row.names=F)
}

finalOut = paste(results.path, 'allResults.Rdata', sep='/')
save(bs_filt, conv_filt, flag_filt, file= finalOut)
