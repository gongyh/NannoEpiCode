#!/bin/bash

#cat HC_VLC.TMM.TPM.matrix | grep -v Genes | sort -k2,2nr | cut -f1,2 > HC.TMM.TPM.matrix
#cat HC_VLC.TMM.TPM.matrix | grep -v Genes | sort -k3,3nr | cut -f1,3 > VLC.TMM.TPM.matrix

#python3.6 splitBed2.py genes2.bed HC.TMM.TPM.matrix HC2
#python3.6 splitBed2.py genes2.bed VLC.TMM.TPM.matrix LC2

'''
computeMatrix reference-point --referencePoint TSS \
        --regionsFileName HC2_G1.bed HC2_G2.bed HC2_G3.bed HC2_G4.bed \
        --scoreFileName H0.mRp.clN.Fnor.smooth.bigWig \
        --outFileName HC.computeMatrix.mat.gz \
        --outFileNameMatrix HC.computeMatrix.vals.mat.gz \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --samplesLabel HC \
        --sortRegions no \
        --numberOfProcessors 32 --binSize 500

computeMatrix reference-point --referencePoint TSS \
        --regionsFileName LC2_G1.bed LC2_G2.bed LC2_G3.bed LC2_G4.bed \
        --scoreFileName H24.mRp.clN.Fnor.smooth.bigWig \
        --outFileName LC.computeMatrix.mat.gz \
        --outFileNameMatrix LC.computeMatrix.vals.mat.gz \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --samplesLabel LC \
        --sortRegions no \
        --numberOfProcessors 32 --binSize 500

plotProfile --matrixFile HC.computeMatrix.mat.gz \
    --outFileName MNase_HC_heatmap.pdf \
    --regionsLabel G1 G2 G3 G4 --plotWidth 10 --plotHeight 9 

plotProfile --matrixFile LC.computeMatrix.mat.gz \
    --outFileName MNase_LC_heatmap.pdf \
    --regionsLabel G1 G2 G3 G4 --plotWidth 10 --plotHeight 9
'''

#python3.6 count_nucleosomes.py > HC_nucleosomes.txt
#python3.6 count_nucleosomes2.py > LC_nucleosomes.txt

python3.6 count_nucleosomes_2k.py > HC_nucleosomes2.txt
python3.6 count_nucleosomes2_2k.py > LC_nucleosomes2.txt

python3.6 count_nucleosomes3.py > HC_nucleosomes_TES.txt
python3.6 count_nucleosomes4.py > LC_nucleosomes_TES.txt

python3.6 splitBed3.py genes2.bed H0_H24_DEG.txt DEG
python3.6 count_nucleosomes_DEG.py > DEG_nucleosomes.txt
python3.6 count_nucleosomes_DEG_TES.py > DEG_nucleosomes2.txt

