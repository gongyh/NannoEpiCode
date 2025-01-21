#!/bin/bash

conda activate gimme
export LD_PRELOAD=$CONDA_PREFIX/lib/libgomp.so

# H3k9ac
#cat H0_H3k9ac_R1_peaks.narrowPeak H0_H3k9ac_R2_peaks.narrowPeak > H0_H3k9ac_peaks.narrowPeak
#gimme motifs H0_H3k9ac_peaks.narrowPeak H0_H3k9ac_motifs -g genome.fa --denovo

#cat H24_H3k9ac_R1_peaks.narrowPeak H24_H3k9ac_R2_peaks.narrowPeak > H24_H3k9ac_peaks.narrowPeak
#gimme motifs H24_H3k9ac_peaks.narrowPeak H24_H3k9ac_motifs -g genome.fa --denovo

#gimme motifs consensus/H3k9ac/deseq2/H0_H3k9acvsH24_H3k9ac/H0_H3k9acvsH24_H3k9ac.deseq2.FDR0.01.results.bed H3k9ac_DEPeaks_motifs -g genome.fa --denovo

# H3k27ac
#cat H0_H3k27ac_R1_peaks.narrowPeak H0_H3k27ac_R2_peaks.narrowPeak > H0_H3k27ac_peaks.narrowPeak
#gimme motifs H0_H3k27ac_peaks.narrowPeak H0_H3k27ac_motifs -g genome.fa --denovo

#cat H24_H3k27ac_R1_peaks.narrowPeak H24_H3k27ac_R2_peaks.narrowPeak > H24_H3k27ac_peaks.narrowPeak
#gimme motifs H24_H3k27ac_peaks.narrowPeak H24_H3k27ac_motifs -g genome.fa --denovo

#gimme motifs consensus/H3k27ac/deseq2/H0_H3k27acvsH24_H3k27ac/H0_H3k27acvsH24_H3k27ac.deseq2.FDR0.01.results.bed H3k27ac_DEPeaks_motifs -g genome.fa --denovo

# Kcr

#cat H0_H3kcr_R1_peaks.narrowPeak H0_H3kcr_R2_peaks.narrowPeak > H0_H3kcr_peaks.narrowPeak
#gimme motifs H0_H3kcr_peaks.narrowPeak H0_H3kcr_motifs -g genome.fa --denovo

#cat H24_H3kcr_R1_peaks.narrowPeak H24_H3kcr_R2_peaks.narrowPeak > H24_H3kcr_peaks.narrowPeak
#gimme motifs H24_H3kcr_peaks.narrowPeak H24_H3kcr_motifs -g genome.fa --denovo

#gimme motifs consensus/H3kcr/deseq2/H0_H3kcrvsH24_H3kcr/H0_H3kcrvsH24_H3kcr.deseq2.FDR0.01.results.bed H3kcr_DEPeaks_motifs -g genome.fa --denovo

# H3k4me2

cat H0_H3k4m2_R1_peaks.narrowPeak H0_H3k4m2_R2_peaks.narrowPeak > H0_H3k4m2_peaks.narrowPeak
gimme motifs H0_H3k4m2_peaks.narrowPeak H0_H3k4m2_motifs -g genome.fa --denovo

cat H24_H3k4m2_R1_peaks.narrowPeak H24_H3k4m2_R2_peaks.narrowPeak > H24_H3k4m2_peaks.narrowPeak
gimme motifs H24_H3k4m2_peaks.narrowPeak H24_H3k4m2_motifs -g genome.fa --denovo

#gimme motifs consensus/H3k4m2/deseq2/H0_H3k4m2vsH24_H3k4m2/H0_H3k4m2vsH24_H3k4m2.deseq2.FDR0.01.results.bed H3k4m2_DEPeaks_motifs -g genome.fa --denovo

