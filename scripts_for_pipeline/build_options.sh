#!/bin/bash
#builds a options file for input into FastOrtho
#usage = build_options.sh project_prefix option_file_name.txt
#eg: build_options.sh run1 run1_option_file.txt

args=("$@")

RESULTS_FILE="${args[0]}.end"


echo "--project_name ${args[0]}
--working_directory $(pwd)
--formatdb_path $MAKEBLASTDB
--blastall_path $BLASTP
--mcl_path $WORK/mcl/bin/mcl
--pv_cutoff $EVALUE
--pi_cutoff 0.0
--pmatch_cutoff 0.0
--maximum_weight 316.0
--result_file ./$RESULTS_FILE
--inflation 1.5
--blast_cpus 12
--blast_b 1000
--blast_v 1000
--blast_e $EVALUE" > ${args[1]}

for faa in $(ls $FAAS/*.fas)
do echo "--single_genome_fasta $faa" >> ${args[1]}
done