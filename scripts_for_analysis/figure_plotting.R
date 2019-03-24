#figure_plotting.R

#LOAD LIBS AND SET GLOBAL VARS
library(tidyverse)
library(ggplot2)
library(cowplot)
setwd("~/gitreps/convergent_evo_coral")
source('scripts_for_analysis/coral_reproduction_functions.R')
xangle=10


##LOAD THE DATASET TO PLOT:
ll=load('overlap_results/moderate/allResults.Rdata') #moderate cutoff set


#Revise clade names in pairings 
conv_filt <- conv_filt %>%
	mutate(pairing = sub('Anti', 'Sis.', cp)) %>%
	mutate(pairing = sub('Anti', 'Sis.', pairing))



#-----------------------------------------------------------------------------------------------#
####1. breakdown of substitition type for verticals by gene

order = conv_filt %>% 
	group_by(pairing) %>%
	summarize(N=n()) %>%
	arrange(desc(N))
order = factor(order$pairing, levels=order$pairing)

vOnly <- conv_filt %>% 
	mutate(pairing=factor(pairing, levels=order)) %>%
	rename(`Change type` = type) %>%
	filter(conv_pair=='VV') %>%
	ggplot() +
		geom_bar(aes(x=pairing, fill=`Change type`)) +
		labs(x='Lineage Pair', y='Overlapping substititions') +
		coord_flip() +
		theme(axis.text.x= element_text(angle=xangle))
		
###1S breakdown for all
all<-conv_filt %>% 
	mutate(pairing=factor(pairing, levels=order)) %>%
	rename(`Change type` = type) %>%
	ggplot() +
		geom_bar(aes(x=pairing, fill=`Change type`)) +
		labs(x='Lineage Pair', y='Overlapping substititions') +
		coord_flip() +
		theme(axis.text.x= element_text(angle= xangle))

plot(all)

#-----------------------------------------------------------------------------------------------#
### S2 BS passing genes + GO (VV and HH)

#set up names 
phylo.ordered.names = c('All vertical', 'All horizontal sisters', 'Montipora', "Montipora's sister", 'Galaxia', "Galaxia's sister", 'Porites', "Porites's sister", 'Pocilloporid', "Pocilloporid's sister")
bsf<-bs_filt %>% 
	mutate(`Significant` = 
		(bs_phenotype=='Vertical' & bs_tpass) | 
		(bs_phenotype=='Horizontal' & bs_apass) ) %>%
	mutate(clade=as.character(clade)) %>%
	mutate(clade = if_else(clade=='AllTarget', 'All vertical', clade)) %>%
	mutate(clade = if_else(clade=='AntiAllTarget', 'All horizontal sisters', clade)) %>%
	mutate(clade = if_else(grepl('Anti', clade), paste(clade, 'sister', sep="'s "), clade)) %>%
	mutate(clade = sub('Anti', '', clade)) %>%
	mutate(clade = factor(clade, levels= rev(phylo.ordered.names)))

#plot counts for the verticals
bsf1 <-	bsf %>%
	filter(bs_phenotype=='Vertical') %>%
	ggplot() + 
	geom_bar(aes(x=clade, fill=`Significant`)) +
	labs(x='Lineage', y='Gene count', subtitle='Vertical') +
	coord_flip() +
	theme(axis.text.x= element_text(angle=10),
		legend.position="none")

#plot counts for horizontals
bsf2 <-	bsf %>%
	filter(bs_phenotype=='Horizontal') %>%
	ggplot() + 
	geom_bar(aes(x=clade, fill=`Significant`)) +
	labs(x='Lineage', y='Gene count', subtitle='Horizontal') +
	coord_flip() +
	theme(axis.text.x= element_text(angle=10),
		legend.position="none")

#gather all counts for each phenotype

tot_genes <- bsf %>%
	group_by(bs_phenotype) %>%
	summarize(tot_genes =length(unique(ortholog)))

tot_any_sig <- bsf %>%
	filter(bs_pass) %>%
	group_by(bs_phenotype) %>%
	summarize(tot_bs_genes =length(unique(ortholog)))

bsf3 <- inner_join(tot_genes, tot_any_sig, by = 'bs_phenotype') %>%
	mutate(not_bs=tot_genes-tot_bs_genes) %>%
	gather(key='pos', value='count', tot_bs_genes, not_bs) %>%
	mutate(`Positive selection`=if_else(pos=='not_bs', FALSE, TRUE)) %>%
	mutate(pheno = if_else(bs_phenotype=='Vertical', 'Vert.', 'Horiz.')) %>%
	ggplot() +
		geom_bar(aes(x=pheno, y=count, fill=`Positive selection`), stat='identity') +
		labs(x='Phenotype', y='Genes')


#plot in grid
# plot_grid(bsf1, bsf2, bsf3, labels=c('A', 'B', 'C'), nrow=1)

#-----------------------------------------------------------------------------------------------#
###S3 convergence frequency comparison
#Do convergence events occur more often in vertical pairs?

#FORMAT PLOTTING DATA

#gather percentages of convergence passing subs
pctdf <-conv_filt %>% 
	mutate(conv_pass = 
		if_else(conv_pass,
			'convergence',
			'no convergence'
			)) %>%
	mutate(conv_pass = factor(conv_pass, levels=c('convergence', 'no convergence'))) %>%
	rename(`Change type`= conv_pass) %>%
	group_by(pairing, conv_pair, `Change type`) %>%
	summarize(N=n()) %>% 
	mutate(pct_conv = N/sum(N)*100) %>%
	filter(`Change type`=='convergence') %>%
	# rename(`Phenotype Pair` = conv_pair) %>%
	arrange(conv_pair, pct_conv)

#get order and colors
ordered = rev(unique(pctdf$pairing))
#set colors
nVV = length(unique(conv_filt$cp[conv_filt$conv_pair=='VV']))
nVH = length(unique(conv_filt$cp[conv_filt$conv_pair=='VH']))
nHH = length(unique(conv_filt$cp[conv_filt$conv_pair=='VH']))
boxplot<-ggplot(pctdf, aes(x=conv_pair, y=pct_conv, fill=conv_pair))
colset<-unique(ggplot_build(boxplot)$data[[1]]['fill'])
axCols = c(rep(colset[3, 'fill'], nVV), rep(colset[2, 'fill'], nVH), rep(colset[1, 'fill'], nHH))


#plot counts by clade pairing
countConv <- conv_filt %>%
	rename(`Convergence`=conv_pass) %>%
	mutate(pairFact = factor(pairing, levels=ordered)) %>%
	ggplot() +
		geom_bar(aes(x= pairFact, fill=`Convergence`)) + 
		labs(x='Lineage pair', y='Overlapping Substitutions') +
		coord_flip() +
		theme(
			axis.text.x= element_text(angle=10),
			axis.text.y = element_text(colour = axCols, face='bold'),
			legend.position=c(0.6,0.8))
		

#plot percentage by clade pairing
pctConv <- pctdf %>%
	mutate(pairFact=factor(pairing, levels=ordered)) %>%
	rename(`Phenotype Pair` = conv_pair) %>%
	ggplot() +
		geom_bar(aes(x= pairFact, y= pct_conv, fill=`Phenotype Pair`), stat='identity') +
		labs(x='', y='Molecular convergence (%)') +
		coord_flip() +
		theme(
			legend.position='none',
			axis.text.y=element_blank())

#plot proportions for convergence filter passing
countConv2<-conv_filt %>% 
	group_by(conv_pair, conv_pass) %>%
	summarize(N=n()) %>%
	mutate(`Percent Convergence`=N/sum(N)*100) %>%
	filter(conv_pass) %>% 
	rename(`Phenotype Pair` = conv_pair) %>%
	ggplot() +
		geom_bar(aes(x=`Phenotype Pair`, y=`Percent Convergence`, fill=`Phenotype Pair`), stat='identity') +
		labs(y='Convergence (% overlapping changes)', x='Phenotype pair') +
		theme(legend.position='none')


#plot mean percentage by phenotype pairing
pctConv2 <- pctdf %>% 
	rename(`Phenotype Pair` = conv_pair) %>%
	rename(`Percent with convergence` = pct_conv) %>%
	ggplot() +
	geom_boxplot(aes(x=`Phenotype Pair`, y=`Percent with convergence`, fill=`Phenotype Pair`)) +
	theme()

#BUILD PLOT
legPos='right'
plot_grid(countConv, pctConv, pctConv2, nrow=1, rel_widths = c(1.5, 1, 1), align = "hv", axis = "tb", labels=LETTERS[1:3])



#STATS FOR THESE PLOTS

#for cumulative counts (barplot)
compareList = list('VH', 'HH')

doFisher = function(x){
	ft_in <- conv_filt %>%
		filter(conv_pair %in% c(x, 'VV')) %>%
		mutate(conv_pair = factor(conv_pair, levels=c('VV', x))) %>%
		mutate(convergent = factor(conv_pass, levels=c(TRUE, FALSE))) %>%
		select(c('conv_pair', 'convergent'))
	fisher.test(ft_in$conv_pair, ft_in$convergent)
}

fish_res = lapply(compareList, doFisher)
names(fish_res) = c('VV_vs_VH', 'VV_vs_HH')
fish_res[['VV_vs_VH']]
fish_res[['VV_vs_HH']]



#for mean by clade pair (boxplot)
#anova
fit <- aov(pct_conv ~ conv_pair, data=pctdf)
summary(fit)

#t.tests
t.test(pct_conv ~ conv_pair, data=pctdf[pctdf$conv_pair %in% c('VV', 'VH'),])
t.test(pct_conv ~ conv_pair, data=pctdf[pctdf$conv_pair %in% c('VV', 'HH'),])




#NUBMERS OF CONVERGENCE EVENTS PER GENE
head(conv_filt)
noDup = conv_filt[!duplicated(conv_filt$site),]
noDup = noDup[noDup$conv_pass & noDup$conv_pair=='VV',]
cSitePerGene <- noDup %>% count(ortholog) %>% as_tibble
allGenes = unique(bs_filt$ortholog[bs_filt$bs_phenotype=='Vertical'])
noOverlapGenes = allGenes[!allGenes %in% cSitePerGene$ortholog]
length(noOverlapGenes)
noSitePerGene = tibble(ortholog=noOverlapGenes, n=0)

mean(cSitePerGene$n)
median(cSitePerGene$n)

#histogram for number of convergence events per gene with any overlapping subs
cSitePerGene %>%
	ggplot() +
		geom_histogram(aes(x=n)) +
		labs(x='Number of convergent sites', y='Count')

#histogram for number of convergence events per gene for all genes
allPerGene = rbind(cSitePerGene, noSitePerGene)
mean(allPerGene$n)
median(allPerGene$n)
hall<-allPerGene %>%
	ggplot() +
		geom_histogram(aes(x=n)) +
		labs(x='Convergence per gene', y='Count')

plot_grid(vOnly+theme(legend.position=c(0.4, 0.77)), hall, rel_widths=c(2,1), labels=LETTERS[1:2])

#-----------------------------------------------------------------------------------------------#
###S3 conv passing genes + GO	
#plot frequency of genes with convergence events

tot_genes <- conv_filt %>%
	group_by(cp, conv_pair) %>%
	summarize(tot=length(unique(ortholog)))

tot_conv_genes <- conv_filt %>%
	filter(conv_pass) %>%
	group_by(cp, conv_pair) %>%
	summarize(conv=length(unique(ortholog)))

gene_df <- inner_join(tot_genes, tot_conv_genes, by='cp')


vv <- gene_df %>%
	mutate(not_conv=tot-conv) %>%
	gather(key='convergenceCall', value='count', conv, not_conv) %>%
	mutate(`Covnergence detected` = if_else(convergenceCall=='conv', TRUE, FALSE)) %>%
	filter(conv_pair.x=='VV') %>%
	arrange(tot)
	
hh <- gene_df %>%
	mutate(not_conv=tot-conv) %>%
	gather(key='convergenceCall', value='count', conv, not_conv) %>%
	mutate(`Covnergence detected` = if_else(convergenceCall=='conv', TRUE, FALSE)) %>%
	filter(conv_pair.x=='HH') %>%
	arrange(tot) %>%
	mutate(pairing = sub('Anti', 'Sister ', cp)) %>%
	mutate(pairing = sub('Anti', 'Sister ', cp))

vv$fcp = factor(vv$cp, levels=rev(unique(vv$cp)))
hh$fcp = factor(hh$pairing, levels=rev(unique(hh$pairing)))

conv1 <- vv %>%
	ggplot() +
		geom_bar(aes(x= fcp, y=count, fill=`Covnergence detected`), stat='identity') +
		labs(y='Genes with overlap', x='Pair', subtitle='Vertical pairs') +
		coord_flip() +
		theme(axis.text.x= element_text(angle=10),
		legend.position="none")
conv2 <- hh %>%
	ggplot() +
		geom_bar(aes(x= fcp, y=count, fill=`Covnergence detected`), stat='identity') +
		labs(y='Genes with overlap', x='Pair', subtitle='Horizontal pairs') +
		coord_flip() +
		theme(axis.text.x= element_text(angle=10),
		legend.position="none")


#set up genes for full set of comparisons
tot_vv = unique(bs_filt$ortholog[bs_filt$bs_phenotype=='Vertical'])
tot_hh = unique(bs_filt$ortholog[bs_filt$bs_phenotype=='Horizontal'])
cvv = unique(conv_filt$ortholog[conv_filt$conv_pair=='VV'])
chh = unique(conv_filt$ortholog[conv_filt$conv_pair=='HH'])
nvv = tot_vv[!tot_vv %in% cvv]
nhh = tot_hh[!tot_hh %in% chh]	

conv3 <- tibble(pairing = factor(c('VV', 'VV', 'HH', 'HH'), levels=c('VV','HH')),
			`Covnergence detected` = c(TRUE, FALSE, TRUE, FALSE),
			count = c(length(cvv), length(nvv), length(chh), length(nhh))) %>%
			ggplot() +
				geom_bar(aes(x=pairing, y=count, fill=`Covnergence detected`), stat='identity') +
				labs(x='Phenotype pairing', y='Gene count')

plot_grid(conv1, conv2, conv3, labels=LETTERS[1:3], nrow=1)


#-----------------------------------------------------------------------------------------------#
###Figure 3 convergence and branch sites genes
##compare frequency of convergence events in positive

#build set of all genes tests



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
clades = unique(bs_filt$clade)
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
	mutate(cp=factor(cp, levels=cp)) %>%
	mutate(pairing = sub('Anti', 'Sis.', cp)) %>%
	mutate(pairing = sub('Anti', 'Sis.', pairing))


pcdat <- mdat %>% 
	mutate(pairing = sub('Anti', 'Sis.', cp)) %>%
	mutate(pairing = sub('Anti', 'Sis.', pairing)) %>%
	mutate(Nnot = nPossible - Ngene) %>%
	select(pairing, Ngene, Nnot, nPossible, conv_pair) %>%
	gather(key=`Conv. & Pos.`, value = 'count', Ngene, Nnot) %>%
	mutate(`Conv. & Pos.` = if_else(`Conv. & Pos.`=='Ngene', TRUE, FALSE)) %>%
	arrange(conv_pair, nPossible) %>%
	mutate(pairFact = factor(pairing, levels=rev(unique(pairing))))


#BUILD PLOTS
#for VV	
pc1 <- pcdat %>%
	filter(conv_pair=='VV') %>%
	ggplot() +
		geom_bar(aes(x= pairFact, y=count, fill=`Conv. & Pos.`, col=`Conv. & Pos.`), stat='identity', size=0.1) +
		guides(size=F) +
		labs(x='Lineage pair', y='Gene count') +
		scale_fill_manual(values=c('gray95', 'grey5')) +
		scale_color_manual(values=c('grey5', 'grey5')) +
		coord_flip() +
		theme(
			axis.text.x= element_text(angle=10),
			# legend.position=c(0.75,0.86),
			legend.position='none',
			plot.margin = unit(c(1,3,1,1), "cm"))
pc1

pc1 + theme(legend.position=c(0.75,0.86))




#for HH
pc2 <- pcdat %>%
	filter(conv_pair=='HH') %>%
	ggplot() +
		geom_bar(aes(x= pairFact, y=count, fill=`Conv. & Pos.`), stat='identity') +
		labs(x='Lineage pair', y='Gene count') +
		coord_flip() +
		theme(
			axis.text.x= element_text(angle=10),
			legend.position=c(0.75,0.86),
			plot.margin = unit(c(1,3,1,1), "cm"))
pc2



#Set up data for cumulative

#get total possible
#(this is all genes with at least two from within a phenotype)
avv<-bs_filt %>%
	filter(bs_phenotype=='Vertical') %>%
	count(ortholog) %>%
	filter(n>2) %>%
	select(ortholog)
ahh<-bs_filt %>%
	filter(bs_phenotype=='Horizontal') %>%
	count(ortholog) %>%
	filter(n>2) %>%
	select(ortholog)

#count unique genes when this occured
cpc <- conv_filt %>% 
	filter(conv_pass) %>%
	filter(c1PosGene | c2PosGene) %>%
	filter(conv_pair %in% c('VV', 'HH')) %>%
	group_by(conv_pair) %>%
	summarize(Nboth=length(unique(ortholog)))
cpc$possible = c(nrow(avv), nrow(ahh))
cpc$Nnot = cpc$possible - cpc$Nboth

pc3 <- cpc %>%
	gather(key=cnp, value='count', Nboth, Nnot) %>%
	mutate(`Conv. & Pos.`=if_else(cnp =='Nboth', TRUE, FALSE)) %>%
	ggplot() +
		geom_bar(aes(x=conv_pair, y=count, fill=`Conv. & Pos.`), stat='identity') +
		labs(x='Phenotype pair', y='Gene count')
	


plot_grid(pc1, pc2, pc3, nrow=1)






#-----------------------------------------------------------------------------------------------#
###Figure Convergence and positive selection substititions comparison

# #set colors
# nVV = length(unique(conv_filt$cp[conv_filt$conv_pair=='VV']))
# nVH = length(unique(conv_filt$cp[conv_filt$conv_pair=='VH']))
# nHH = length(unique(conv_filt$cp[conv_filt$conv_pair=='VH']))
# colset<-ggplot_build(pctConv2)$data[[1]]['fill']
# axCols = c(rep(colset[3, 'fill'], nVV), rep(colset[2, 'fill'], nVH), rep(colset[1, 'fill'], nHH))

order =unique(mdat$pairing)

#get the proportion of convergent substitions that occured in positive genes vs not positive genes
pc_filt <- conv_filt %>%
	mutate(both = conv_pass & (c1PosGene | c2PosGene)) %>%
	group_by(pairing, conv_pair, both) %>%
	summarize(count = n()) %>%
	mutate(pct = count/sum(count)*100) %>%
	arrange(conv_pair, pct)


order = rev(unique(pc_filt$pairing[pc_filt$both]))
pc_filt<-pc_filt %>%
	rename(`Conv. & Pos.`=both) %>%
	mutate(pairFact = factor(pairing, levels=order))

pcCompare1 <- pc_filt %>%
	ggplot() +
		geom_bar(aes(x= pairFact, y=count, fill=`Conv. & Pos.`), stat='identity') +
		labs(x='Lineage pair', y='Overlapping substition count') +
		coord_flip() +
		theme(axis.text.y = element_text(colour = axCols, face='bold'), 
			axis.text.x = element_text(angle= xangle),
			legend.position='none')


#plot percentage by convergence pairing
pcCompare2<-pc_filt %>%
	filter(`Conv. & Pos.`) %>%
	rename(`Phenotype pair`=conv_pair) %>%
	ggplot() +
		geom_bar(aes(x= pairFact, y=pct, fill=`Phenotype pair`), stat='identity') +
		labs(x='', y='Conv. & Pos. (%)') +
		coord_flip() +
		theme(legend.position='none',
			  axis.text.y = element_blank())

#plot percentage boxplots
pcCompare3<-pc_filt %>%
	filter(`Conv. & Pos.`) %>%
	rename(`Phenotype pair`=conv_pair) %>%
	ggplot() +
		geom_boxplot(aes(x=`Phenotype pair`, y=pct, fill=`Phenotype pair`)) +
		labs(x='Phenotype pair', y='Conv. & Pos. (%)')
 box_legend = get_legend(pcCompare3)



plot_grid(pcCompare1+theme(legend.position=c(0.6,0.8)), pcCompare2,  pcCompare3, nrow=1, rel_widths = c(1.5, 1, 1), align = "hv", axis = "tb", labels=c(LETTERS[1:3]))

#-----------------------------------------------------------------------------------------------#
###Figure 4 compare frequency of convergence and flagged substititions

pc_cf <- conv_filt %>%
	mutate(both = conv_pass & (c1FlagSite | c2FlagSite)) %>%
	group_by(pairing, conv_pair, both) %>%
	summarize(count = n()) %>%
	mutate(pct = count/sum(count)*100) %>%
	arrange(conv_pair, both, pairing) %>%
	mutate(pctBoth = 100-pct) %>%
	arrange(conv_pair, pctBoth) %>%
	as_tibble

o=pc_cf %>% filter(!both) %>% arrange(conv_pair, pctBoth)
order = rev(unique(o$pairing))


cfCompare1 <- pc_cf %>%
	rename(`Conv. & Pos.Site`=both) %>%
	mutate(pairFact = factor(pairing, levels=order)) %>%
	ggplot() +
		geom_bar(aes(x= pairFact, y=count, fill=`Conv. & Pos.Site`), stat='identity') +
		labs(x='Lineage pair', y='Overlapping substition count') +
		coord_flip() +
		theme(axis.text.y = element_text(colour = axCols, face='bold'), 
			axis.text.x = element_text(angle= xangle),
			legend.position=c(0.475,0.7))


#plot percentage by convergence pairing
cfCompare2<-pc_cf %>%
	filter(!both) %>%
	mutate(pairFact = factor(pairing, levels=order)) %>%
	rename(`Phenotype pair`=conv_pair) %>% 
	ggplot() +
		geom_bar(aes(x= pairFact, y=pctBoth, fill=`Phenotype pair`), stat='identity') +
		labs(x='', y='Conv. & Pos. sites (%)') +
		coord_flip() +
		theme(
			legend.position='none',
			axis.text.y=element_blank())
		
vvhh = pctdf %>% filter(conv_pair != 'VH')
	t.test(vvhh$pct_conv~vvhh$conv_pair)


#plot percentage boxplots
cfCompare3<-pc_cf %>%
	filter(!both) %>%
	rename(`Phenotype pair`=conv_pair) %>%
	ggplot() +
		geom_boxplot(aes(x=`Phenotype pair`, y= pctBoth, fill=`Phenotype pair`)) +
		labs(x='Phenotype pair', y='Conv. & Pos. site (%)')	


plot_grid(cfCompare1, cfCompare2, cfCompare3, rel_widths =c(1.6, 1, 1), align = "hv", axis = "tb", labels=c(LETTERS[1:2]), nrow=1)





#------------------------------------------ frequency comparisons for each cutoff set









