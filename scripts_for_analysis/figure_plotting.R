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



# Figure 2 breakdown of substitition type for verticals by gene ---------

#get order for barplot
order = conv_filt %>% 
	group_by(pairing) %>%
	summarize(N=n()) %>%
	arrange(desc(N))
order = factor(order$pairing, levels=order$pairing)

#build barplot (2A)
f2a <- conv_filt %>% 
	mutate(pairing=factor(pairing, levels=order)) %>%
	rename(`Change type` = type) %>%
	filter(conv_pair=='VV') %>%
	ggplot() +
		geom_bar(aes(x=pairing, fill=`Change type`)) +
		labs(x='Lineage Pair', y='Overlapping substititions') +
		coord_flip() +
		theme(axis.text.x= element_text(angle=xangle),
		      legend.position=c(.4, 0.78)) +
    guides(fill=guide_legend(title="Change type (example)"))

#get number of convergence events per gene
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


#build histogram for number of convergence events per gene for genes with at least one convergence event
#(Fig 2B)
convGeneHist=cSitePerGene %>%
  ggplot() +
  geom_histogram(aes(x=n)) +
  labs(x='Convergent sites per gene', y='Count')


#build same histogram for all genes (Figure 2B)
allPerGene = rbind(cSitePerGene, noSitePerGene)
f2b<-allPerGene %>%
  ggplot() +
  geom_histogram(aes(x=n)) +
  labs(x='Convergence per gene', y='Count')

mean(allPerGene$n)
median(allPerGene$n)


plot_grid(f2a, f2b)


# Fig S2 breakdown of change type for all ---------------------------------

all<-conv_filt %>% 
	mutate(pairing=factor(pairing, levels=order)) %>%
	rename("Change type" = type) %>%
	ggplot() +
		geom_bar(aes(x=pairing, fill=`Change type`)) +
		labs(x='Lineage Pair', y='Overlapping substititions') +
		coord_flip() +
		theme(axis.text.x= element_text(angle= xangle),
		      legend.position=c(.375, 0.7)) +
    guides(fill=guide_legend(title="Change type")) +
    scale_y_continuous(limits=c(0,41000))

plot(all)


# Figure 4 convergence frequency comparison ------------------------------
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
  guides(fill=guide_legend(title="Transmission\nmode\npair")) +
  labs(x='Transmission mode pair', y='Percent convergence') 

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


# Figure 3 convergence among positive selected genes ----------------------
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


#BUILD PLOT
#for VV	
f3a <- pcdat %>%
	filter(conv_pair=='VV') %>%
	ggplot() +
		geom_bar(aes(x= pairFact, y=count, fill=`Conv. & Pos.`, col=`Conv. & Pos.`), stat='identity', size=0.1, color='black') +
		labs(x='Lineage pair', y='Gene count') +
		scale_fill_manual(values=c('gray95', 'grey5')) +
		# scale_color_manual(values=c('grey5', 'grey5')) +
		coord_flip() +
		theme(
			axis.text.x= element_text(angle=10),
			plot.margin = unit(c(1,3,1,1), "cm"),
			legend.position=c(0.75,0.86)) +
    guides(fill=guide_legend(title="Conv. & Pos."))
f3a


# figure S5  ----------------------------------------------------
#Equivalent to figure 3 but done for H-H pairs
#for HH
pc2 <- pcdat %>%
	filter(conv_pair=='HH') %>%
	ggplot() +
		geom_bar(aes(x= pairFact, y=count, fill=`Conv. & Pos.`), stat='identity', size=0.1, color='black') +
		labs(x='Lineage pair', y='Gene count') +
		coord_flip() +
    scale_fill_manual(values=c('gray95', 'grey5')) +
		theme(
			axis.text.x= element_text(angle=10),
			legend.position=c(0.75,0.86),
			plot.margin = unit(c(1,3,1,1), "cm")) +
    guides(fill=guide_legend(title="Conv. & Pos."))
pc2



# Figure S3 Convergence and positive selection substititions comparis --------

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
			legend.position=c(0.575,0.775)) +
  guides(fill=guide_legend(title="Conv. & Pos."))


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
		labs(x='Transmission mode pairing', y='Conv. & Pos. (%)') +
    guides(fill=guide_legend(title="Transmission\nmode\npairing"))

box_legend = get_legend(pcCompare3)



plot_grid(pcCompare1, pcCompare2,  pcCompare3, nrow=1, rel_widths = c(1.5, 1, 1), align = "hv", axis = "tb", labels=c(LETTERS[1:3]))



# Figure S4 compare frequency of convergence and flagged substitutions --------

pc_cf <- conv_filt %>%
	mutate(both = conv_pass & (c1FlagSite | c2FlagSite) ) %>%
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
			legend.position=c(0.5,0.75)) +
    guides(fill=guide_legend(title="Conv. & Pos.Site"))


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
		labs(x='Transmission mode pairing', y='Conv. & Pos. site (%)')	+
    guides(fill=guide_legend(title="Transmission\nmode\npairing"))

quartz()
plot_grid(cfCompare1, cfCompare2, cfCompare3, rel_widths =c(1.6, 1, 1), align = "hv", axis = "tb", labels=c(LETTERS[1:2]), nrow=1)



#double-check the VV convergent and flagged genes
flaggedConvGenes = conv_filt %>%
  mutate(both = conv_pass & (c1FlagSite | c2FlagSite) ) %>% 
  filter(both,
         conv_pair=='VV') %>% 
  pull(ortholog) %>% 
  unique() %>% 
  length()
