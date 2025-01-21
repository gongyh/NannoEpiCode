#!/bin/bash

#ln -s /mnt/data7/gongyh/Nanno/epi/histone/IMET1v2/Sequence/WholeGenomeFasta/genome.fa genome.fa
#ln -s /mnt/data7/gongyh/Nanno/epi/histone/IMET1v2/Annotation/Genes/genes.gtf genes.gtf

annotatePeaks.pl rna genome.fa -p H0_H3k4m2_R1_summits.bed -gid -gtf genes.gtf -cpu 12 > H0_H3k4m2_R1_summits.annotatePeaks.txt
annotatePeaks.pl rna genome.fa -p H0_H3k4m2_R2_summits.bed -gid -gtf genes.gtf -cpu 12 > H0_H3k4m2_R2_summits.annotatePeaks.txt
annotatePeaks.pl rna genome.fa -p H24_H3k4m2_R1_summits.bed -gid -gtf genes.gtf -cpu 12 > H24_H3k4m2_R1_summits.annotatePeaks.txt
annotatePeaks.pl rna genome.fa -p H24_H3k4m2_R2_summits.bed -gid -gtf genes.gtf -cpu 12 > H24_H3k4m2_R2_summits.annotatePeaks.txt

