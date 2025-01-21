#!/bin/bash

names=("H0_H3k9ac_R1" "H24_H3k9ac_R1" "H0_H3k27ac_R1" "H24_H3k27ac_R1" "H0_H3k4m2_R1" "H24_H3k4m2_R1" "H0_H3kcr_R1" "H24_H3kcr_R1" "H0_H3k9ac_R2" "H24_H3k9ac_R2" "H0_H3k27ac_R2" "H24_H3k27ac_R2" "H0_H3k4m2_R2" "H24_H3k4m2_R2" "H0_H3kcr_R2" "H24_H3kcr_R2")

for name in ${names[*]}; do
  plotProfile --matrixFile ../${name}.computeMatrix.mat.gz \
        --outFileName ${name}.plotProfile.pdf \
        --outFileNameData ${name}.plotProfile.tab
'''
  plotHeatmap --matrixFile ${name}.computeMatrix.mat.gz --kmeans 5 \
        --outFileName ${name}_k5.plotHeatmap.pdf --heatmapWidth 6 \
        --outFileNameMatrix ${name}_k5.plotHeatmap.mat.tab \
        --outFileSortedRegions ${name}_k5.plotHeatmap.sorted.bed \
        --zMin auto auto --zMax auto auto
'''
done

