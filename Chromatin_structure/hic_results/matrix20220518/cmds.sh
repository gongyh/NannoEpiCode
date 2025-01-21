#!/bin/bash

## assign genes to bins
bedtools intersect -a bins_10000.bed -b IMET1v2.geneOnly.bed -wo -F 0.5 | cut -f1,2,3,7 > bins_10000_genes.bed
bedtools intersect -a bins_100000.bed -b IMET1v2.geneOnly.bed -wo -F 0.5 | cut -f1,2,3,7 > bins_100000_genes.bed
bedtools intersect -a bins_40000.bed -b IMET1v2.geneOnly.bed -wo -F 0.5 | cut -f1,2,3,7 > bins_40000_genes.bed

awk -F'\t' 'BEGIN{print "chrom_start_end_GD"}{seg=$1"_"$2"_"$3;gd[seg]+=1}END{for(s in gd){print s"_"gd[s]}}' bins_40000_genes.bed | sed 's/_/\t/g' | awk 'NR == 1; NR > 1 {print $0 | "sort -k1,1V -k2,2n"}' /dev/stdin > bins_40000_gd.txt

cat H3k9ac.consensus_peaks.results.txt | awk -F'\t' 'BEGIN{OFS="\t"}NR>1{print $2,$3,$4,$1}' > H3k9ac.consensus_peaks.bed
cat H3k27ac.consensus_peaks.results.txt | awk -F'\t' 'BEGIN{OFS="\t"}NR>1{print $2,$3,$4,$1}' > H3k27ac.consensus_peaks.bed
cat H3kcr.consensus_peaks.results.txt | awk -F'\t' 'BEGIN{OFS="\t"}NR>1{print $2,$3,$4,$1}' > H3kcr.consensus_peaks.bed
cat H3k4m2.consensus_peaks.results.txt | awk -F'\t' 'BEGIN{OFS="\t"}NR>1{print $2,$3,$4,$1}' > H3k4m2.consensus_peaks.bed

bedtools intersect -a bins_40000.bed -b H3k9ac.consensus_peaks.bed -wo -F 0.5 | cut -f1,2,3,7 > bins_40000_H3k9ac.bed
bedtools intersect -a bins_40000.bed -b H3k27ac.consensus_peaks.bed -wo -F 0.5 | cut -f1,2,3,7 > bins_40000_H3k27ac.bed
bedtools intersect -a bins_40000.bed -b H3kcr.consensus_peaks.bed -wo -F 0.5 | cut -f1,2,3,7 > bins_40000_H3kcr.bed
bedtools intersect -a bins_40000.bed -b H3k4m2.consensus_peaks.bed -wo -F 0.5 | cut -f1,2,3,7 > bins_40000_H3k4m2.bed

