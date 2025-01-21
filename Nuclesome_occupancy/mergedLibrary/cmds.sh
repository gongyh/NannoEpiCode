#!/bin/bash

#bamCoverage -b H0_R1.mLb.clN.sorted.bam -o H0_R1.mLb.clN.bw --MNase --binSize 1 -p 36
#bamCoverage -b H0_R2.mLb.clN.sorted.bam -o H0_R2.mLb.clN.bw --MNase --binSize 1 -p 36
#bamCoverage -b H24_R1.mLb.clN.sorted.bam -o H24_R1.mLb.clN.bw --MNase --binSize 1 -p 36
#bamCoverage -b H24_R2.mLb.clN.sorted.bam -o H24_R2.mLb.clN.bw --MNase --binSize 1 -p 36

bamCoverage -b H0_R1.mLb.clN.sorted.bam -o H0_R1.mLb.clN.bw --binSize 1 -p 36 
bamCoverage -b H0_R2.mLb.clN.sorted.bam -o H0_R2.mLb.clN.bw --binSize 1 -p 36
bamCoverage -b H24_R1.mLb.clN.sorted.bam -o H24_R1.mLb.clN.bw --binSize 1 -p 36
bamCoverage -b H24_R2.mLb.clN.sorted.bam -o H24_R2.mLb.clN.bw --binSize 1 -p 36

bedtools makewindows -b genome.bed -w 500 | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,"bin_"NR}'  > genome_500.bed

bigWigAverageOverBed H0_R1.mLb.clN.bw genome_500.bed H0_R1.mLb.clN.tab
bigWigAverageOverBed H0_R2.mLb.clN.bw genome_500.bed H0_R2.mLb.clN.tab
bigWigAverageOverBed H24_R1.mLb.clN.bw genome_500.bed H24_R1.mLb.clN.tab
bigWigAverageOverBed H24_R2.mLb.clN.bw genome_500.bed H24_R2.mLb.clN.tab

Rscript pearsonCmp_bin500.R

