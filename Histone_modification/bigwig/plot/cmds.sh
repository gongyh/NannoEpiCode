#!/bin/bash

conda activate pygenometracks

pyGenomeTracks --tracks tracks_poor.ini --region chr1:1350000-1550000 --width 20 --outFileName gene_poor.pdf
pyGenomeTracks --tracks tracks_CAH1.ini --region chr20:209500-214500 --width 20 --outFileName gene_CAH1.pdf

pyGenomeTracks --tracks tracks_ME1.ini --region chr26:168200-174300 --width 15 --outFileName gene_ME1.pdf --fontSize 6

pyGenomeTracks --tracks tracks_RNAseq.ini --BED genome.fa.include_regions.bed --width 20 --outFileName RNAseq.pdf

bigwigCompare -bs 1 -b1 H24.Fnor.smooth.bigWig -b2 H0.Fnor.smooth.bigWig --pseudocount 0.15 -p 16 -of bigwig -o log2fc.Fnor.smooth.bigWig
bigwigCompare -bs 1 -b1 VLC24.bw -b2 0h.bw --pseudocount 0.06 -p 16 -of bigwig -o log2fc.RNAseq.bw

pyGenomeTracks --tracks tracks_PMA2.ini --region chr17:162000-168300 --width 15 --outFileName NoPMA2_HCLC.pdf
