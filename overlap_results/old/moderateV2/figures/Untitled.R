
bs_ctrim = bs_filt %>% 
	mutate(clade = as.character(clade)) %>%
	select(clade, ortholog, bs_pass) 
	
conv_bs <- conv_filt %>% 
   left_join(
     bs_ctrim %>% rename(b1Clade = clade, b1_bs_passable = bs_pass), 
     by = c("b1Clade", "ortholog")) %>%
   left_join(
     bs_ctrim %>% rename(b2Clade = clade, b2_bs_passable = bs_pass), 
         by = c("b2Clade", "ortholog")) %>% 
   as_tibble
     

sum(conv_bs$ortholog==conv_filt$ortholog)==nrow(conv_filt)
sum(conv_bs$b1_bs_passable ==conv_filt$c1PosGene)==nrow(conv_filt)
conv_filt[conv_filt$c1PosGene != conv_bs$b1_bs_passable,]
conv_bs[conv_filt$c1PosGene != conv_bs$b1_bs_passable,] 
     
     
     
bs_ftrim <- bs_filt %>% 
	mutate(bs_passable = p.values < bs.gene.pvalue.cutoff) %>%
	select(clade, ortholog, bs_passable)

flag_filt0 <- flag_dat %>%
	left_join(bs_ftrim, by = c('clade', 'ortholog')) %>%
	as_tibble

flag_filt0