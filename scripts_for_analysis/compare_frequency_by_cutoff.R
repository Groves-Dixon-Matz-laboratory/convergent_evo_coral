

#LOAD LIBS AND SET GLOBAL VARS
library(tidyverse)
library(ggplot2)
library(cowplot)
setwd("~/gitreps/convergent_evo_coral")
source('scripts_for_analysis/coral_reproduction_functions.R')
xangle=10

##LOAD THE DATASET TO PLOT:
# ll=load('overlap_results/relaxed/allResults.Rdata') #relaxed cutoff set
# ll=load('overlap_results/relaxed_no_brood/allResults.Rdata') #moderate cutoff set
# ll=load('overlap_results/moderate/allResults.Rdata') #moderate cutoff set
# ll=load('overlap_results/moderate_no_brood/allResults.Rdata') #moderate no brood cutoff set
# ll=load('overlap_results/semi_strict/allResults.Rdata') #semistrict no brood cutoff set
# ll=load('overlap_results/strict/allResults.Rdata') #semistrict no brood cutoff set

cutoff_sets = c('relaxed', 'relaxed_no_brood', 'moderate', 'moderate_no_brood', 'semi_strict', 'strict')

#RUN INCLUDING VH

for (co in cutoff_sets){
	print('----------------')
	print(paste('Running', co))
	fileName = paste('overlap_results', co, 'allResults.Rdata', sep='/')
	load(fileName)


	#Revise clade names in pairings 
	conv_filt <- conv_filt %>%
		mutate(pairing = sub('Anti', 'Sis.', cp)) %>%
		mutate(pairing = sub('Anti', 'Sis.', pairing))
	
	
	
	#-------- convergence frequency
	
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
	
	#plot mean percentage by phenotype pairing
	#basic stats
	fit <- aov(pct_conv ~ conv_pair, data=pctdf)
	attributes(summary(fit))
	pval=round(summary(fit)[[1]][5][1,1], digits=3)
	
	
	convBox <- pctdf %>% 
		rename(`Phenotype Pair` = conv_pair) %>%
		rename(`Percent convergence` = pct_conv) %>%
		ggplot() +
			geom_boxplot(aes(x=`Phenotype Pair`, y=`Percent convergence`, fill=`Phenotype Pair`)) +
			labs(title=co, subtitle=paste("aov p-value =", pval)) +
			theme(legend.position='none')
	
	
	#-------- convergence in branch site gene frequency
	
	#get the proportion of convergent substitions that occured in positive genes vs not positive genes
	pc_filt <- conv_filt %>%
		mutate(both = conv_pass & (c1PosGene | c2PosGene)) %>%
		group_by(pairing, conv_pair, both) %>%
		summarize(count = n()) %>%
		mutate(pct = count/sum(count)*100) %>%
		arrange(conv_pair, pct) %>%
		filter(both)
	
	#plot percentage boxplots
	#basic stats
	fit <- aov(pct ~ conv_pair, data= pc_filt)
	pval=round(summary(fit)[[1]][5][1,1], digits=3)
	
	
	pcBox<-pc_filt %>%
		filter(both) %>%
		rename(`Phenotype pair`=conv_pair) %>%
		ggplot() +
			geom_boxplot(aes(x=`Phenotype pair`, y=pct, fill=`Phenotype pair`)) +
			labs(title=co, x='Phenotype pair', y='Conv. & Pos. (%)', subtitle=paste("aov p-value =", pval)) +
			theme(legend.position='none')
	
	
	#-------- convergence and flagged
	
	
	pc_cf <- conv_filt %>%
		mutate(both = conv_pass & (c1FlagSite | c2FlagSite)) %>%
		group_by(pairing, conv_pair, both) %>%
		summarize(count = n()) %>%
		mutate(pct = count/sum(count)*100) %>%
		arrange(conv_pair, both, pairing) %>%
		mutate(pctBoth = 100-pct) %>%
		arrange(conv_pair, pctBoth) %>%
		filter(!both) %>%
		as_tibble
	
	
	#plot percentage boxplots
	#basic stats
	fit <- aov(pctBoth ~ conv_pair, data= pc_cf)
	summary(fit)
	pval=round(summary(fit)[[1]][5][1,1], digits=3)
	
	
	#plot percentage boxplots
	fcBox<-pc_cf %>%
		filter(!both) %>%
		rename(`Phenotype pair`=conv_pair) %>%
		ggplot() +
			geom_boxplot(aes(x=`Phenotype pair`, y= pctBoth, fill=`Phenotype pair`)) +
			labs(title=co, x='Phenotype pair', y='Conv. & Pos. site (%)', subtitle=paste("aov p-value =", pval))	
			
			
	plot(plot_grid(convBox, pcBox, fcBox, nrow=1, rel_widths=c(1,1,1.4)))
}






#RUN NOT INCLUDING VH

for (co in cutoff_sets){
	print('----------------')
	print(paste('Running', co))
	fileName = paste('overlap_results', co, 'allResults.Rdata', sep='/')
	load(fileName)
	conv_filt = conv_filt[conv_filt$conv_pair != 'VH',]
	conv_filt$conv_pair = factor(conv_filt$conv_pair, levels=c('VV', 'HH'))


	#Revise clade names in pairings 
	conv_filt <- conv_filt %>%
		mutate(pairing = sub('Anti', 'Sis.', cp)) %>%
		mutate(pairing = sub('Anti', 'Sis.', pairing))
	
	
	
	#-------- convergence frequency
	
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
	
	#plot mean percentage by phenotype pairing
	#basic stats
	ttest = t.test(pctdf$pct_conv~ pctdf$conv_pair)
	pval=round(ttest$p.value, digits=3)
	
	
	convBox <- pctdf %>% 
		rename(`Phenotype Pair` = conv_pair) %>%
		rename(`Percent convergence` = pct_conv) %>%
		ggplot() +
			geom_boxplot(aes(x=`Phenotype Pair`, y=`Percent convergence`, fill=`Phenotype Pair`)) +
			labs(title=co, subtitle=paste("t-test p-value =", pval)) +
			theme(legend.position='none')
	
	
	#-------- convergence in branch site gene frequency
	
	#get the proportion of convergent substitions that occured in positive genes vs not positive genes
	pc_filt <- conv_filt %>%
		mutate(both = conv_pass & (c1PosGene | c2PosGene)) %>%
		group_by(pairing, conv_pair, both) %>%
		summarize(count = n()) %>%
		mutate(pct = count/sum(count)*100) %>%
		arrange(conv_pair, pct) %>%
		filter(both)
	
	#plot percentage boxplots
	#basic stats
	ttest = t.test(pc_filt$pct~ pc_filt$conv_pair)
	pval=round(ttest$p.value, digits=3)
	
	
	pcBox<-pc_filt %>%
		filter(both) %>%
		rename(`Phenotype pair`=conv_pair) %>%
		ggplot() +
			geom_boxplot(aes(x=`Phenotype pair`, y=pct, fill=`Phenotype pair`)) +
			labs(title=co, x='Phenotype pair', y='Conv. & Pos. (%)', subtitle=paste("t-test p-value =", pval)) +
			theme(legend.position='none')
	
	
	#-------- convergence and flagged
	
	
	pc_cf <- conv_filt %>%
		mutate(both = conv_pass & (c1FlagSite | c2FlagSite)) %>%
		group_by(pairing, conv_pair, both) %>%
		summarize(count = n()) %>%
		mutate(pct = count/sum(count)*100) %>%
		arrange(conv_pair, both, pairing) %>%
		mutate(pctBoth = 100-pct) %>%
		arrange(conv_pair, pctBoth) %>%
		filter(!both) %>%
		as_tibble
	
	
	#plot percentage boxplots
	#basic stats
	ttest = t.test(pc_cf$pctBoth ~ pc_cf$conv_pair)
	pval=round(ttest$p.value, digits=3)
	
	
	#plot percentage boxplots
	fcBox<-pc_cf %>%
		filter(!both) %>%
		rename(`Phenotype pair`=conv_pair) %>%
		ggplot() +
			geom_boxplot(aes(x=`Phenotype pair`, y= pctBoth, fill=`Phenotype pair`)) +
			labs(title=co, x='Phenotype pair', y='Conv. & Pos. site (%)', subtitle=paste("t-test p-value =", pval))	
			
			
	plot(plot_grid(convBox, pcBox, fcBox, nrow=1, rel_widths=c(1,1,1.4)))
}

