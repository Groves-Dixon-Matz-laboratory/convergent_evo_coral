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


# basic model -------------------------------------------------------------
basic_model = glm(genConvergent~conv_pair, family='binomial', data=conv_filt)
summary(basic_model)



# basic vv and vh only ----------------------------------------------------
basic_vh_model = glm(genConvergent~conv_pair, family='binomial', data=conv_vh)
summary(basic_vh_model)


# random effect of gene ---------------------------------------------------
null_rgene_model = glmer(genConvergent~(1|ortholog), family='binomial', data=conv_filt)
alt_rgene_model = glmer(genConvergent~conv_pair + (1|ortholog), family='binomial', data=conv_filt)
summary(null_rgene_model)
summary(alt_rgene_model)
anova(null_rgene_model, alt_rgene_model, test="Chisq")

# random effect of gene VV and VH only---------------------------------------------------
null_rgene_vh_model = glmer(genConvergent~(1|ortholog), family='binomial', data=conv_vh)
alt_rgene_vh_model = glmer(genConvergent~conv_pair + (1|ortholog), family='binomial', data=conv_vh)
summary(null_rgene_vh_model)
summary(alt_rgene_vh_model)
anova(null_rgene_vh_model, alt_rgene_vh_model)



# random effect for pairing -----------------------------------------------

null_rpair_model = glmer(genConvergent~ (1|pairing), family='binomial', data=conv_filt)
alt_rpair_model = glmer(genConvergent~conv_pair + (1|pairing), family='binomial', data=conv_filt)
summary(null_rpair_model)
summary(alt_rpair_model)
anova(null_rpair_model, alt_rpair_model)



# format for lineage effects ----------------------------------------------
#get all lineages
lins = unique(append(conv_filt$b1Clade, conv_filt$b2Clade))

#set up new df with binaries for each lineage
ldat = conv_filt %>% 
  select(-b1post, -b2post, -type, -b1Anc, -b2Anc, -gcp, -cp, -b1New, -b2New)
for (l in lins){
  print(paste(l, '..', sep='.'))
  ldat[[l]] = (ldat$b1Clade==l) | (ldat$b2Clade==l)
}

#create a vh only equivalent
vhldat = ldat %>% 
  filter(conv_pair !='HH') %>% 
  mutate(conv_pair=factor(conv_pair, levels=c('VH', 'VV')))


#check correlations for dummy assignments
jldat = ldat[,15:ncol(ldat)]
cor(jldat)


# fixed effect for lineages------------------------------------------
null_fixedLin_model = glm(genConvergent ~ 
            AntiPocilloporid +
            Pocilloporid +
            Galaxia +
            Montipora +
            AntiMontipora +
            AntiGalaxia +
            AntiPorites +
            Porites,
          family='binomial', data=ldat)

alt_fixedLin_model = glm(genConvergent ~ conv_pair +
            AntiPocilloporid +
            Pocilloporid +
            Galaxia +
            Montipora +
            AntiMontipora +
            AntiGalaxia +
            AntiPorites +
            Porites,
          family='binomial', data=ldat)
summary(null_fixedLin_model)
alias(null_fixedLin_model)
summary(alt_fixedLin_model)
anova(null_fixedLin_model, alt_fixedLin_model, test="Chisq")


# repeat with singularities removed ---------------------------------------
null_fixedLin_model = glm(genConvergent ~ 
                            Pocilloporid +
                            Galaxia +
                            Montipora +
                            Porites,
                          family='binomial', data=ldat)

alt_fixedLin_model = glm(genConvergent ~ conv_pair +
                           Pocilloporid +
                           Galaxia +
                           Montipora +
                           Porites,
                         family='binomial', data=ldat)
summary(null_fixedLin_model)
alias(null_fixedLin_model)
summary(alt_fixedLin_model)
anova(null_fixedLin_model, alt_fixedLin_model, test="Chisq")



# fixed interaction effect for lineages------------------------------------------
fixedIntLin_model = glm(genConvergent ~ conv_pair +
                          conv_pair*Pocilloporid +
                          conv_pair*AntiPocilloporid +
                          conv_pair*Galaxia +
                          conv_pair*AntiGalaxia +
                          conv_pair*Montipora +
                          conv_pair*AntiMontipora +
                          conv_pair*Porites +
                          conv_pair*AntiPorites,
                         family='binomial', data=vhldat)
summary(fixedIntLin_model)
alias(fixedIntLin_model)

# repeat fixed interaction effect for lineages remove singluarities------------------------------------------
fixedIntLin_model = glm(genConvergent ~ conv_pair +
                          conv_pair*Pocilloporid +
                          conv_pair*Galaxia +
                          conv_pair*Montipora,
                        family='binomial', data=vhldat)
summary(fixedIntLin_model)


# interactive random effect for lineage -----------------------------------------------
randIntLin_model = glmer(genConvergent ~ conv_pair +
               (1+conv_pair|AntiPocilloporid) +
               (1+conv_pair|Pocilloporid) +
               (1+conv_pair|Galaxia) +
               (1+conv_pair|Montipora) +
               (1+conv_pair|AntiMontipora) +
               (1+conv_pair|AntiGalaxia) +
               (1+conv_pair|AntiPorites) +
               (1+conv_pair|Porites),
             family='binomial', data=ldat)
summary(randIntLin_model)


# interactive random effect for lineage vh only -----------------------------------------------
randIntLin_vh_model = glmer(genConvergent ~ conv_pair +
                           (1+conv_pair|AntiPocilloporid) +
                           (1+conv_pair|Pocilloporid) +
                           (1+conv_pair|Galaxia) +
                           (1+conv_pair|Montipora) +
                           (1+conv_pair|AntiMontipora) +
                           (1+conv_pair|AntiGalaxia) +
                           (1+conv_pair|AntiPorites) +
                           (1+conv_pair|Porites),
                         family='binomial', data=vhldat)
summary(randIntLin_vh_model)

# interaction model -------------------------------------------------------
intMod = glm(genConvergent ~ conv_pair +
            conv_pair:AntiPocilloporid +
            conv_pair:Pocilloporid +
            conv_pair:Galaxia +
            conv_pair:Montipora +
            conv_pair:AntiMontipora +
            conv_pair:AntiGalaxia +
            conv_pair:AntiPorites +
            conv_pair:Porites,
          family='binomial', data=ldat)
summary(intMod)
