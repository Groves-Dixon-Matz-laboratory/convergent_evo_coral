#convergence_glms.R
#use glms to test significance for variation in frequecy of convergence events between phyenotype pairs
#hypothesis is that V-V pairs (both vertical transmitters) have higher rate of convergence events per
#overlapping substititutions (indepdent amino acid changes at same position)


#libs and wd
library(tidyverse)
library(ggplot2)
library(cowplot)
library(lme4)
setwd("~/gitreps/convergent_evo_coral")
source('scripts_for_analysis/coral_reproduction_functions.R')


#upload dataset
co='moderate'
print(paste('Running', co))
fileName = paste('overlap_results', co, 'allResults.Rdata', sep='/')
load(fileName)


#Revise clade names in pairings 
conv_filt <- conv_filt %>%
  mutate(pairing = sub('Anti', 'Sis.', cp)) %>%
  mutate(pairing = sub('Anti', 'Sis.', pairing)) %>% 
  mutate(conv_pair=factor(conv_pair, levels=c('HH', 'VH', 'VV')))

#Set up dataset with only VV and VH
conv_vh = conv_filt %>% 
  filter(conv_pair != 'HH') %>% 
  mutate(conv_pair=factor(conv_pair, levels=c('VH', 'VV')))


#---------- GLMS
#first stack the pairs to get single lineage per line
sub1 = conv_filt %>% 
  select(ortholog, site, b1Clade, conv_pair, genConvergent) %>% 
  rename('clade'=b1Clade)
sub2 = conv_filt %>% 
  select(ortholog, site, b2Clade, conv_pair, genConvergent) %>% 
  rename('clade'=b2Clade)
sldat = rbind(sub1, sub2)
#simple model with gene as random effect
library(lme4)
sldat$conv_pair = factor(sldat$conv_pair, levels=c('HH', 'VH', 'VV'))
null_basic_model = glmer(genConvergent~(1|ortholog), family='binomial', data=sldat)
alt_basic_model = glmer(genConvergent~conv_pair + (1|ortholog), family='binomial', data=sldat)
summary(null_basic_model)
summary(alt_basic_model)
anova(null_basic_model, alt_basic_model)


#models for frequency of convergence events
base_mod = glm(genConvergent~1,
               family='binomial',
               data=sldat)

gene_only = glmer(genConvergent~
                    (1|ortholog),
                  family='binomial',
                  data=sldat)

lin_only = glmer(genConvergent~
                   (1|clade),
                 family='binomial',
                 data=sldat)

lin_gene = glmer(genConvergent~
                   (1|ortholog) +
                   (1|clade),
                 family='binomial',
                 data=sldat)

all_three = glmer(genConvergent~
                    conv_pair +
                    (1|ortholog) +
                    (1|clade),
                  family='binomial',
                  data=sldat)

interaction = glmer(genConvergent~
                      conv_pair +
                      (1|ortholog) +
                      (1+conv_pair|clade),
                    family='binomial',
                    data=sldat)

#model summaries
summary(gene_only)
summary(lin_only)
summary(lin_gene)
summary(all_three)
summary(interaction)

#compare
anova(lin_only, lin_gene)
anova(gene_only, lin_gene)
anova(gene_only, all_three)
anova(lin_gene, all_three)
anova(all_three, interaction)

ll=load('~/Desktop/mods.Rdata')
save(lin_only, gene_only, lin_gene, all_three, interaction, file='~/Desktop/mods.Rdata')