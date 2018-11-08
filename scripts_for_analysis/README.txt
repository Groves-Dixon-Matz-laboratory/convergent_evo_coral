#scripts_for_analysis
#This directory contains scripts used for downstream analysis after 
#the ortholog pipeline and PAML have been run. The main input files
#for these can be found in ancestral_reconstruction/ and branch_sites_tests/.


#Contents:
coral_reproduction_functions.R  -- functions used by other R scripts
LRT_for_branch_sites.R          -- likelihood ratio test on null and alternative likelihoods from branch site tests
prepareConvergenceInput.R       -- take in the raw ancestral reconstruction data (output from parse_multi_rstV2.py) and prepare for input into branch_sites_overlap.R
overlap_analysis.R              -- assesses levels of overlap between Branch site tests, individual sites flagged in the BS tests, and convergence events
	moderate_cutoff_set.R           -- set of moderate cutoffs used in overlap_analysis.R
	relaxed_cutoff_set.R            -- set of relaxed cutoffs used in overlap_analysis.R
figure_plotting.R	            -- plots final figures for publication
summarizing.R                   -- used to get numbers for manuscript such as total number of sites showing convergence
compare_frequency_by_cutoff.R   -- build boxplots comparing frequency of convergence between phenotypes for each cutoff set