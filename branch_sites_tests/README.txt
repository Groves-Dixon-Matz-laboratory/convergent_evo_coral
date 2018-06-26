#branch_sites_tests
#This directory contains the files for the branch sites tests
#10 Branch sites tests were performed:
#one for each vertical transmitting clade individually (4 total)
#one for each horizontal transmitting clade individually (4 total)
#one for the full set of vertical trasmitters
#one for the full set of horizontal transmitters


#contents:
#null likelihood files:
	nullLikelihoods_branchSites_*.tsv
		These have the likelihoods for the null model for each gene
#alternative likelihood files:
	altLikelihoods_branchSites_*.tsv
		These have the likelihoods for the alternative model.
		Compare the null and alternative likelihoods with LRT_for_branch_sites.R
#likelihood ratio test results:
	bs_lrt*.tsv
		These have the p-values from the likelihood ratio tests
#significant site files:
	*sigBsSites.tsv
		These files have all the sites flagged as significant in the CODEML output files
