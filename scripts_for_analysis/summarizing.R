library(tidyverse)
setwd("~/gitreps/convergent_evo_coral")

ll= load('~/gitreps/convergent_evo_coral/overlap_results/moderate/allResults.Rdata')

#total convergence events between verticals
pass = conv_filt %>% filter(conv_pair=='VV' & conv_pass) 
length(unique(pass$site))
length(unique(pass$ortholog))

#total convergence events between verticals
conv = conv_filt %>% filter(conv_pair=='VV' & conv_pass & type=='convergent (L>V; S>V)') 
para = conv_filt %>% filter(conv_pair=='VV' & conv_pass & type=='parallel (L>V; L>V)') 
tc=length(unique(conv$site))
tp=length(unique(para$site))
tc
tp
tc+tp



#number of genes showing at least one VV convergence event
pass = conv_filt %>% filter(conv_pair=='VV' & conv_pass)
length(unique(pass$ortholog))



#total branch-site passing
head(bs_filt)
gids = read.table("ortholog_tables/singleCopyAnnotationsWithDescriptions.tsv", sep='\t', quote="", header = T, stringsAsFactors=F) %>% as_tibble
colnames(gids)[1]='ortholog'
passingGenes = unique(bs_filt$ortholog[bs_filt$bs_pass & bs_filt$bs_phenotype=='Vertical'])
length(passingGenes)
bsPassTable<- bs_filt %>% 
	filter(bs_pass, bs_phenotype=='Vertical') %>%
	left_join(gids, by = 'ortholog') %>%
	mutate(number=as.numeric(sub('ORTHOMCL', '', ortholog))) %>%
	arrange(number) %>%
	select(ortholog, p.values, adj.p, clade, swissProt, egn, db, description)
colnames(bsPassTable) = c('ortholog', 'raw_pvalue', 'adjusted_p', 'clade', 'swissprot', 'entrez', 'entrez_db', 'description')
write.table(bsPassTable, file='~/Desktop/branch_site_passing_genes.tsv', sep='\t', quote=F, row.names=F)	










pass = conv_filt %>% mutate(posPass=c1PosGene | c2PosGene) %>% filter(conv_pair=='VV' & conv_pass & posPass)
length(unique(pass$ortholog))



conv_filt %>% filter(ortholog=='ORTHOMCL5637' & conv_pass) %>% arrange(pos) 
conv_filt %>% filter(ortholog=='ORTHOMCL5637' & conv_pass) %>% arrange(pos) %>% select(site) %>% unique()

conv_filt %>% filter(ortholog=='ORTHOMCL2985' & conv_pass) %>% arrange(pos)