#!/bin/bash

#conda activate pygenometracks
#bigwigCompare -bs 1 -b1 H24.Fnor.smooth.bigWig -b2 H0.Fnor.smooth.bigWig --pseudocount 0.15 -p 16 -of bigwig -o log2fc.Fnor.smooth.bigWig
#bigwigCompare -bs 1 -b1 VLC24.bw -b2 0h.bw --pseudocount 0.06 -p 16 -of bigwig -o log2fc.RNAseq.bw

#pyGenomeTracks --tracks tracks_NO18G02680.ini --region chr18:868937-876568 --width 20 --outFileName NO18G02680.pdf
pyGenomeTracks --tracks tracks_NO18G02680.ini --region chr18:868937-876568 --width 20 --outFileName NO18G02680_hq.pdf

pyGenomeTracks --tracks tracks_fig1.ini --region chr18:868937-876568 --width 15 --height 20 --outFileName mnase_fig1.pdf

