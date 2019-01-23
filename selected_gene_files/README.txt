selected_gene_files

#this directory contains intermediate pipeline files for the three genes discussed in manuscript:

ABL1 (ORTHOMCL8234)
filamin C (ORTHOMCL8658)
poly(rC) binding protein 2 (ORTHOMCL8545)


#branch-site test files are separated into subdirectories labeled by vertically transmitting clades

#files types (prefixed by ortholog name):
*.aln                 -- protein alignment
*.codon               -- reverse translated codon sequences
*.codon.tree          -- gene tree constructed from codon sequences
*.rst                 -- ancestral reconstruction output file
*.rstrodTree.newick   -- tree from ancestral reconstruction output
*_ANCESTRAL_RECON.cnt -- ancestral reconstruction control file
*_NULL.cnt            -- control file for null model for branch-site test
*_ALT.cnt             -- control file for alternative model for branch-site test
*_NULL.codeml         -- output file for branch-site null model
*_ALT.codeml          -- output for branch-site alternative model
*_CDS.fasta           -- fasta file of codings sequences from Transdecoder
*_PRO.fasta           -- fasta of protein sequences from Transdecoder (input for MAFFT)
