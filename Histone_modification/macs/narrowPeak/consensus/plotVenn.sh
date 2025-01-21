#!/bin/bash

'''
(.venv)[gongyh@gnode7 consensus]$ wc -l */*.bed
   6130 H3k27ac/H3k27ac.consensus_peaks.bed
   9987 H3k4m2/H3k4m2.consensus_peaks.bed
   6228 H3k9ac/H3k9ac.consensus_peaks.bed
   6202 H3kcr/H3kcr.consensus_peaks.bed
   28547 total
'''

'''
bedtools intersect -wa -a H3k9ac/H3k9ac.consensus_peaks.bed -b H3k27ac/H3k27ac.consensus_peaks.bed > H3k9ac_H3k27ac.bed
bedtools intersect -wa -a H3k9ac/H3k9ac.consensus_peaks.bed -b H3kcr/H3kcr.consensus_peaks.bed > H3k9ac_H3kcr.bed

bedtools intersect -wa -a H3k27ac/H3k27ac.consensus_peaks.bed -b H3k9ac/H3k9ac.consensus_peaks.bed > H3k27ac_H3k9ac.bed
bedtools intersect -wa -a H3k27ac/H3k27ac.consensus_peaks.bed -b H3kcr/H3kcr.consensus_peaks.bed > H3k27ac_H3kcr.bed

bedtools intersect -wa -a H3kcr/H3kcr.consensus_peaks.bed -b H3k9ac/H3k9ac.consensus_peaks.bed > H3kcr_H3k9ac.bed
bedtools intersect -wa -a H3kcr/H3kcr.consensus_peaks.bed -b H3k27ac/H3k27ac.consensus_peaks.bed > H3kcr_H3k27ac.bed
'''

'''
(.venv)[gongyh@gnode7 consensus]$ wc -l *.bed
   6055 H3k27ac_H3k9ac.bed
      6107 H3k27ac_H3kcr.bed
         6055 H3k9ac_H3k27ac.bed
    6125 H3k9ac_H3kcr.bed
       6107 H3kcr_H3k27ac.bed
          6125 H3kcr_H3k9ac.bed
    36574 total
    (.venv)[gongyh@gnode7 consensus]$ cat H3k9ac_H3k27ac.bed H3k9ac_H3kcr.bed | sort | uniq -d | wc -l
    6033
    (.venv)[gongyh@gnode7 consensus]$ cat H3k27ac_H3k9ac.bed H3k27ac_H3kcr.bed | sort | uniq -d | wc -l
    5898
    (.venv)[gongyh@gnode7 consensus]$ cat H3kcr_H3k9ac.bed H3kcr_H3k27ac.bed | sort | uniq -d | wc -l
    5939
'''

cat H3k9ac/H3k9ac.consensus_peaks.boolean.annotatePeaks.txt | cut -f40 | sort | uniq | grep -v "Entrez ID" > H3k9ac_genes.txt
cat H3k27ac/H3k27ac.consensus_peaks.boolean.annotatePeaks.txt | cut -f40 | sort | uniq | grep -v "Entrez ID" > H3k27ac_genes.txt
cat H3kcr/H3kcr.consensus_peaks.boolean.annotatePeaks.txt | cut -f40 | sort | uniq | grep -v "Entrez ID" > H3kcr_genes.txt
cat H3k4m2/H3k4m2.consensus_peaks.boolean.annotatePeaks.txt | cut -f40 | sort | uniq | grep -v "Entrez ID" > H3k4m2_genes.txt


