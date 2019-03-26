#---- CUTOFF SELECTION ----#
set.label = 'moderate'
results.path = 'overlap_results/moderate'


#Pairs to skip for whatever readson
#convergence events between these will be ignored
#(order they are given in the pairs doesn't matter)
#(leave empty if want to keep all)
remove_conv_clade_pairs = list()
remove_conv_clade_spp_pairs = list()


#Branch-sites test filters
bs.ptype = "adj.p"             #type of p value to use (adj.p or pvalues)
bs.alpha = 0.1                 #significance cutoff for chosen p value
bs.min.lineage = 1             #number of target lineages that must be significant to keep gene

#Convergence filters
ancestral.file = "ancestral_reconstruction/ancestral_recon_full_clades_lineages.Rdata" #which ancestral dataset to use. Generate these with prepareConvergenceInput.R
anc.sub.types = c('parallel', 'convergent') #the type of overlapping substitutions that count as convergent
anc.post.cutoff = 0.8                       #the posterior probability cutoff to use to consider ancestral calls
min.convergent.changes = 2                  #the number of lineages that must show a convergence event at a site for it to be kept. Note these don't all have to be the SAME convergent change. For example two clades could show P>V at a position, and two others could show P>I. This would return 4 lineages showing convergence.


#branch-sites flagged sites filters
bs.site.post.cutoff = 0.8        #cutoff for posterior probability the site has a dN/dS > 1
bs.gene.pvalue.cutoff = 0.05     #unadjusted p-value the gene has to have for full branch site test to consider a flagged site passing
min.flagged.tests = 1            #the site must have been flagged in >= this number of tests within a phenotype for the site to pass
