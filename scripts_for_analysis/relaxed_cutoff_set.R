#---- CUTOFF SELECTION ----#
set.label = 'relaxed'
results.path = 'overlap_results/relaxed'

#Branch-sites test filters
bs.ptype = "p.values"             #type of p value to use (adj.p or pvalues)
bs.alpha = 0.05                 #significance cutoff for chosen p value
bs.min.lineage = 1             #number of target lineages that must be significant to keep gene

#Convergence filters
ancestral.file = "ancestral_reconstruction/ancestral_recon_full_clades_lineages.Rdata" #which ancestral dataset to use. Generate these with prepareConvergenceInput.R
anc.sub.types = c('parallel', 'convergent')  #the type of overlapping substitutions that count as convergent
anc.post.cutoff = 0.8                        #the posterior probability cutoff to use to consider ancestral calls
min.convergent = 1                           #the number of convergent pairs that must be seen to count


#branch-sites flagged sites filters
bs.site.post.cutoff = 0.8                    #cutoff for posterior probability the site has a dN/dS > 1



