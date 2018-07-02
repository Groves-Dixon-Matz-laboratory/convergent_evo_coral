#scripts_for_analysis
#This directory contains scripts used for downstream analysis after 
#the ortholog pipeline and PAML have been run. The main input files
#for these can be found in ancestral_reconstruction/ and branch_sites_tests/.


#Contents:
LRT_for_branch_sites.R          -- likelihood ratio test on null and alternative likelihoods from branch site tests
branch_sites_overlap.R          -- subset genes based on overlap between 1. branch site overall significance 2. molecular convergence 3. amino acid position significance from branch site tests
coral_reproduction_functions.R  -- functions used by other R scripts
moderate_cutoff_set.R           -- set of cutoffs used in branch_sites_overlap.R
prepareConvergenceInput.R       -- take in the raw ancestral reconstruction data (output from parse_multi_rstV2.py) and prepare for input into branch_sites_overlap.R
relaxed_cutoff_set.R            -- relaxed set of cutoffs used in branch_sites_overlap.R