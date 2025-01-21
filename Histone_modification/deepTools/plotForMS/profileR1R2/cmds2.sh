#!/bin/bash

'''
clusters=("1" "2" "3" "4" "5")
for cls in ${clusters[*]}; do
  cat ../H0_H3k9ac_k5.plotHeatmap.sorted.bed | grep cluster_${cls} > H0_H3k9ac_c${cls}.bed
  cat ../H0_H3k27ac_k5.plotHeatmap.sorted.bed | grep cluster_${cls} > H0_H3k27ac_c${cls}.bed
  cat ../H0_H3kcr_k5.plotHeatmap.sorted.bed | grep cluster_${cls} > H0_H3kcr_c${cls}.bed
  cat ../H0_H3k4m2_k5.plotHeatmap.sorted.bed | grep cluster_${cls} > H0_H3k4m2_c${cls}.bed
done

names=("H0_H3k9ac" "H24_H3k9ac" "H0_H3k27ac" "H24_H3k27ac" "H0_H3k4m2" "H24_H3k4m2" "H0_H3kcr" "H24_H3kcr")
for name in ${names[*]}; do

    computeMatrix scale-regions \
        --regionsFileName ${name/H24/H0}_c1.bed ${name/H24/H0}_c2.bed ${name/H24/H0}_c3.bed ${name/H24/H0}_c4.bed ${name/H24/H0}_c5.bed \
        --scoreFileName ../${name}_R1.bigWig \
        --outFileName ${name}_R1.computeMatrix.mat.gz \
        --outFileNameMatrix ${name}_R1.computeMatrix.vals.mat.tab \
        --regionBodyLength 5000 \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --smartLabels \
        --numberOfProcessors 48

done
'''

computeMatrixOperations cbind -m H0_H3k9ac_R1.computeMatrix.mat.gz H24_H3k9ac_R1.computeMatrix.mat.gz \
  -o R1_H3k9ac.computeMatrix.mat.gz

computeMatrixOperations cbind -m H0_H3k27ac_R1.computeMatrix.mat.gz H24_H3k27ac_R1.computeMatrix.mat.gz \
  -o R1_H3k27ac.computeMatrix.mat.gz

computeMatrixOperations cbind -m H0_H3k4m2_R1.computeMatrix.mat.gz H24_H3k4m2_R1.computeMatrix.mat.gz \
  -o R1_H3k4m2.computeMatrix.mat.gz

computeMatrixOperations cbind -m H0_H3kcr_R1.computeMatrix.mat.gz H24_H3kcr_R1.computeMatrix.mat.gz \
  -o R1_H3kcr.computeMatrix.mat.gz

names=("R1_H3k9ac" "R1_H3k27ac" "R1_H3k4m2" "R1_H3kcr")

for name in ${names[*]}; do

plotHeatmap --matrixFile ${name}.computeMatrix.mat.gz \
    --outFileName ${name}_k5.plotHeatmap.pdf --heatmapWidth 6 \
    --outFileNameMatrix ${name}_k5.plotHeatmap.mat.tab \
    --zMin auto auto --zMax auto auto --regionsLabel C1 C2 C3 C4 C5

done

