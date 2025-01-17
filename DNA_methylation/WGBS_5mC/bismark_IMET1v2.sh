#!/bin/bash

SNs=("0_1" "0_2" "0_3" "24_1" "24_2" "24_3")

~/progs/Bismark_v0.19.0/bismark_genome_preparation --bowtie2 IMET1v2_bt2 # only need to call once

for SN in ${SNs[*]}; do

  ~/progs/Bismark_v0.19.0/bismark --parallel 48 -o map2 --prefix IMET1v2 --bowtie2 IMET1v2_bt2 \
    -1 Trim/${SN}_R1_val_1.fq.gz -2 Trim/${SN}_R2_val_2.fq.gz

  cd map2

  ~/progs/Bismark_v0.19.0/deduplicate_bismark -p --bam IMET1v2.${SN}_R1_val_1_bismark_bt2_pe.bam

  ~/progs/Bismark_v0.19.0/bismark_methylation_extractor -p --no_overlap --ignore 6 --ignore_r2 6 --comprehensive \
    --report --parallel 48 --bedGraph --zero_based --cutoff 2 --remove_spaces --buffer_size 20G --cytosine_report \
    --CX --genome_folder ../IMET1v2_bt2 IMET1v2.${SN}_R1_val_1_bismark_bt2_pe.deduplicated.bam

  ~/progs/Bismark_v0.19.0/bam2nuc --genome_folder ../IMET1v2_bt2 IMET1v2.${SN}_R1_val_1_bismark_bt2_pe.deduplicated.bam

  #~/progs/Bismark_v0.19.0/bam2nuc --genome_folder ../IMET1v2_bt2 IMET1v2.${SN}_R1_val_1_bismark_bt2_pe.bam

  ~/progs/Bismark_v0.19.0/bismark2report --alignment_report IMET1v2.${SN}_R1_val_1_bismark_bt2_PE_report.txt

  ~/progs/Bismark_v0.19.0/bismark2summary # call for comparison of all samples

  ~/src/Bismark-0.19.0/coverage2cytosine --gzip --genome_folder ../IMET1v2_bt2 \
    -o IMET1v2_${SN}_bismark_bt2_pe.deduplicated.bismark.cov.gz \
    IMET1v2.${SN}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz

  cd ..

done

