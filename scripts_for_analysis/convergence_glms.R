#convergence_glms.R
#use glms to test significance for variation in frequecy of convergence events between phyenotype pairs
#hypothesis is that V-V pairs (both vertical transmitters) have higher rate of convergence events per
#overlapping substititutions (indepdent amino acid changes at same position)


#libs and wd
library(tidyverse)
library(ggplot2)
library(cowplot)
library(lme4)
library(lmerTest)# to get p-value estimations that are not part of the standard lme4 packages
library(sjstats)
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
  mutate(conv_pair=factor(conv_pair, levels=c('HH', 'VH', 'VV')),
         pheno1 = if_else(grepl('Anti', b1Clade),
                          'H',
                          'V'),
         pheno2 = if_else(grepl('Anti', b2Clade),
                          'H',
                          'V')
         )


# basic test for conv pair on convergence -----------------------------------------------
intercept = glm(genConvergent~1,
                family='binomial',
                data=conv_filt)

basic = glm(genConvergent~conv_pair,
               family='binomial',
               data=conv_filt)

#results
llnull = -logLik(intercept)
llalt = -logLik(basic)
x2=2*(llnull - llalt)
summary(basic)
anova(intercept, basic, test='Chisq')
x2  #double-check x2 matches as with glmer

# basic test for convenience_pair on convergence & positive sel -----------
conv_filt = conv_filt %>% 
  mutate(convPos = genConvergent & (c1PosGene | c2PosGene))

interceptCP = glm(convPos~1,
                family='binomial',
                data=conv_filt)

basicCP = glm(convPos~conv_pair,
            family='binomial',
            data=conv_filt)

#results
llnull = -logLik(interceptCP)
llalt = -logLik(basicCP)
x2=2*(llnull - llalt)
summary(basicCP)
anova(interceptCP, basicCP, test='Chisq')
x2


# control for lineage and gene-----------------------------------------------------


#first stack the pairs to get single lineage per line
sub1 = conv_filt %>% 
  mutate(cladeHorizontal = grepl('Anti', b1Clade),
         pairHorizontal = grepl('Anti', b2Clade)) %>% 
  rename('clade'=b1Clade,
         'posGene'=c1PosGene,
         'flagged'=c1FlagSite) %>% 
  select(ortholog, site, clade, cladeHorizontal, pairHorizontal, conv_pair, genConvergent, posGene, flagged) 
sub2 = conv_filt %>% 
  mutate(cladeHorizontal=grepl('Anti', b2Clade),
         pairHorizontal = grepl('Anti', b1Clade)) %>% 
  rename('clade'=b2Clade,
         'posGene'=c2PosGene,
         'flagged'=c2FlagSite) %>% 
  select(ortholog, site, clade, cladeHorizontal, pairHorizontal, conv_pair, genConvergent, posGene, flagged) 
sldat = rbind(sub1, sub2) %>% 
  mutate(convPos = genConvergent & posGene,
         convFlag = genConvergent & flagged,
         vertical = !grepl('Anti', clade),
         conv_pair = factor(conv_pair))
sldat$conv_pair = factor(sldat$conv_pair, levels=c('HH', 'VH', 'VV'))



# models for frequency of convergence events ------------------------------

library(lme4)
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

save(lin_only, gene_only, lin_gene, all_three, interaction, file='~/Desktop/convergence_freq_models.Rdata')
ll=load('~/Desktop/convergence_freq_models.Rdata')


#model summaries
summary(gene_only)
summary(lin_only)
summary(tmode_only)
summary(lin_gene)
summary(all_three)

#compare
anova(lin_only, lin_gene, test='LRT')
anova(gene_only, lin_gene, test='LRT')
anova(gene_only, all_three, test='LRT')
anova(lin_gene, all_three, test='LRT')
llnull = -logLik(lin_gene)
llalt = -logLik(all_three)
x2=2*(llnull - llalt)
x2



# model for co-occurance of convergence and pos selection -----------------


gene_only = glmer(convPos~
                    (1|ortholog),
                  family='binomial',
                  data=sldat)

lin_only = glmer(convPos~
                   (1|clade),
                 family='binomial',
                 data=sldat)

lin_gene = glmer(convPos~
                   (1|ortholog) +
                   (1|clade),
                 family='binomial',
                 data=sldat)

all_three = glmer(convPos~
                    conv_pair +
                    (1|ortholog) +
                    (1|clade),
                  family='binomial',
                  data=sldat)

save(lin_only, gene_only, lin_gene, all_three, interaction, file='~/Desktop/convergence_freq_models.Rdata')
ll=load('~/Desktop/convergence_freq_models.Rdata')


#model summaries
summary(gene_only)
summary(lin_only)
summary(tmode_only)
summary(lin_gene)
summary(all_three)

#compare
anova(lin_only, lin_gene, test='LRT')
anova(gene_only, lin_gene, test='LRT')
anova(gene_only, all_three, test='LRT')
anova(lin_gene, all_three, test='LRT')
llnull = -logLik(lin_gene)
llalt = -logLik(all_three)
x2=2*(llnull - llalt)
x2
