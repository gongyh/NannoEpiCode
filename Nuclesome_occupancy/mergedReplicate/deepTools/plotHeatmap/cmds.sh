#!/bin/bash

'''
cat HC_VLC.TMM.TPM.matrix | grep -v Genes | sort -k2,2nr | cut -f1,2 > HC.TMM.TPM.matrix
cat HC_VLC.TMM.TPM.matrix | grep -v Genes | sort -k3,3nr | cut -f1,3 > VLC.TMM.TPM.matrix

python3.6 splitBed.py genes2.bed HC.TMM.TPM.matrix HC
python3.6 splitBed.py genes2.bed VLC.TMM.TPM.matrix VLC
'''

computeMatrix reference-point --referencePoint TSS \
        --regionsFileName HC_G1.bed HC_G2.bed HC_G3.bed HC_G4.bed \
        --scoreFileName H0.mRp.clN.Fnor.smooth.bigWig \
        --outFileName H0.computeMatrix.mat.gz \
        --outFileNameMatrix H0.computeMatrix.vals.mat.gz \
        --beforeRegionStartLength 500 \
        --afterRegionStartLength 1000 \
        --samplesLabel H0 \
        --sortRegions no \
        --numberOfProcessors 32

computeMatrix reference-point --referencePoint TSS \
        --regionsFileName VLC_G1.bed VLC_G2.bed VLC_G3.bed VLC_G4.bed \
        --scoreFileName H24.mRp.clN.Fnor.smooth.bigWig \
        --outFileName H24.computeMatrix.mat.gz \
        --outFileNameMatrix H24.computeMatrix.vals.mat.gz \
        --beforeRegionStartLength 500 \
        --afterRegionStartLength 1000 \
        --samplesLabel H24 \
        --sortRegions no \
        --numberOfProcessors 32

plotHeatmap --matrixFile H0.computeMatrix.mat.gz \
    --outFileName MNase_H0_heatmap.pdf \
    --sortRegions no --heatmapWidth 6

plotHeatmap --matrixFile H24.computeMatrix.mat.gz \
    --outFileName MNase_H24_heatmap.pdf \
    --sortRegions no --heatmapWidth 6

