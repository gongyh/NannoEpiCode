#!/bin/bash

'''
bigwigCompare -b1 ../../bigwig/H0_H3k9ac_R1.bigWig -b2 ../../bigwig/Input_H0_R1.bigWig -bs 10 -p 48 -o H0_H3k9ac_R1.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H0_H3k27ac_R1.bigWig -b2 ../../bigwig/Input_H0_R1.bigWig -bs 10 -p 48 -o H0_H3k27ac_R1.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H0_H3k4m2_R1.bigWig -b2 ../../bigwig/Input_H0_R1.bigWig -bs 10 -p 48 -o H0_H3k4m2_R1.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H0_H3kcr_R1.bigWig -b2 ../../bigwig/Input_H0_R1.bigWig -bs 10 -p 48 -o H0_H3kcr_R1.bigWig -of bigwig

bigwigCompare -b1 ../../bigwig/H0_H3k9ac_R2.bigWig -b2 ../../bigwig/Input_H0_R2.bigWig -bs 10 -p 48 -o H0_H3k9ac_R2.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H0_H3k27ac_R2.bigWig -b2 ../../bigwig/Input_H0_R2.bigWig -bs 10 -p 48 -o H0_H3k27ac_R2.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H0_H3k4m2_R2.bigWig -b2 ../../bigwig/Input_H0_R2.bigWig -bs 10 -p 48 -o H0_H3k4m2_R2.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H0_H3kcr_R2.bigWig -b2 ../../bigwig/Input_H0_R2.bigWig -bs 10 -p 48 -o H0_H3kcr_R2.bigWig -of bigwig

bigwigCompare -b1 ../../bigwig/H24_H3k9ac_R1.bigWig -b2 ../../bigwig/Input_H24_R1.bigWig -bs 10 -p 48 -o H24_H3k9ac_R1.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H24_H3k27ac_R1.bigWig -b2 ../../bigwig/Input_H24_R1.bigWig -bs 10 -p 48 -o H24_H3k27ac_R1.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H24_H3k4m2_R1.bigWig -b2 ../../bigwig/Input_H24_R1.bigWig -bs 10 -p 48 -o H24_H3k4m2_R1.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H24_H3kcr_R1.bigWig -b2 ../../bigwig/Input_H24_R1.bigWig -bs 10 -p 48 -o H24_H3kcr_R1.bigWig -of bigwig

bigwigCompare -b1 ../../bigwig/H24_H3k9ac_R2.bigWig -b2 ../../bigwig/Input_H24_R2.bigWig -bs 10 -p 48 -o H24_H3k9ac_R2.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H24_H3k27ac_R2.bigWig -b2 ../../bigwig/Input_H24_R2.bigWig -bs 10 -p 48 -o H24_H3k27ac_R2.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H24_H3k4m2_R2.bigWig -b2 ../../bigwig/Input_H24_R2.bigWig -bs 10 -p 48 -o H24_H3k4m2_R2.bigWig -of bigwig
bigwigCompare -b1 ../../bigwig/H24_H3kcr_R2.bigWig -b2 ../../bigwig/Input_H24_R2.bigWig -bs 10 -p 48 -o H24_H3kcr_R2.bigWig -of bigwig

cat genes.bed | grep "\.1" > genes2.bed

names=("H0_H3k9ac" "H24_H3k9ac" "H0_H3k27ac" "H24_H3k27ac" "H0_H3k4m2" "H24_H3k4m2" "H0_H3kcr" "H24_H3kcr")
for name in ${names[*]}; do

computeMatrix scale-regions \
        --regionsFileName genes2.bed \
        --scoreFileName ${name}_R1.bigWig \
        --outFileName ${name}_R1.computeMatrix.mat.gz \
        --outFileNameMatrix ${name}_R1.computeMatrix.vals.mat.tab \
        --regionBodyLength 5000 \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --smartLabels \
        --numberOfProcessors 48

computeMatrix scale-regions \
        --regionsFileName genes2.bed \
        --scoreFileName ${name}_R2.bigWig \
        --outFileName ${name}_R2.computeMatrix.mat.gz \
        --outFileNameMatrix ${name}_R2.computeMatrix.vals.mat.tab \
        --regionBodyLength 5000 \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --smartLabels \
        --numberOfProcessors 48

computeMatrix scale-regions \
        --regionsFileName genes2.bed \
        --scoreFileName ${name}_R1.bigWig ${name}_R2.bigWig \
        --outFileName ${name}.computeMatrix.mat.gz \
        --outFileNameMatrix ${name}.computeMatrix.vals.mat.tab \
        --regionBodyLength 5000 \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --smartLabels \
        --numberOfProcessors 48

done


computeMatrixOperations rbind -m H0_H3k9ac_R1.computeMatrix.mat.gz H0_H3k9ac_R2.computeMatrix.mat.gz \
  -o H0_H3k9ac.computeMatrix.mat.gz
computeMatrixOperations rbind -m H24_H3k9ac_R1.computeMatrix.mat.gz H24_H3k9ac_R2.computeMatrix.mat.gz \
  -o H24_H3k9ac.computeMatrix.mat.gz

computeMatrixOperations rbind -m H0_H3k27ac_R1.computeMatrix.mat.gz H0_H3k27ac_R2.computeMatrix.mat.gz \
  -o H0_H3k27ac.computeMatrix.mat.gz
computeMatrixOperations rbind -m H24_H3k27ac_R1.computeMatrix.mat.gz H24_H3k27ac_R2.computeMatrix.mat.gz \
  -o H24_H3k27ac.computeMatrix.mat.gz

computeMatrixOperations rbind -m H0_H3k4m2_R1.computeMatrix.mat.gz H0_H3k4m2_R2.computeMatrix.mat.gz \
  -o H0_H3k4m2.computeMatrix.mat.gz
computeMatrixOperations rbind -m H24_H3k4m2_R1.computeMatrix.mat.gz H24_H3k4m2_R2.computeMatrix.mat.gz \
  -o H24_H3k4m2.computeMatrix.mat.gz

computeMatrixOperations rbind -m H0_H3kcr_R1.computeMatrix.mat.gz H0_H3kcr_R2.computeMatrix.mat.gz \
  -o H0_H3kcr.computeMatrix.mat.gz
computeMatrixOperations rbind -m H24_H3kcr_R1.computeMatrix.mat.gz H24_H3kcr_R2.computeMatrix.mat.gz \
  -o H24_H3kcr.computeMatrix.mat.gz
'''

names=("H0_H3k9ac" "H24_H3k9ac" "H0_H3k27ac" "H24_H3k27ac" "H0_H3k4m2" "H24_H3k4m2" "H0_H3kcr" "H24_H3kcr")
for name in ${names[*]}; do
'''
  plotProfile --matrixFile ${name}.computeMatrix.mat.gz \
        --outFileName ${name}.plotProfile.pdf \
        --outFileNameData ${name}.plotProfile.tab
'''
  plotHeatmap --matrixFile ${name}.computeMatrix.mat.gz --kmeans 5 \
        --outFileName ${name}_k5.plotHeatmap.pdf --heatmapWidth 6 \
        --outFileNameMatrix ${name}_k5.plotHeatmap.mat.tab \
        --outFileSortedRegions ${name}_k5.plotHeatmap.sorted.bed \
        --zMin auto auto --zMax auto auto
done

