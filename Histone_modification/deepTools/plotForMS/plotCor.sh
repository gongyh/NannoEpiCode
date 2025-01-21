#!/bin/bash

#multiBigwigSummary bins -b ../../bigwig/*.bigWig -o all.npz -bs 200 --smartLabels -p 48

#plotCorrelation -in all.npz -c pearson -p scatterplot -o cor_scatter.pdf --skipZeros --removeOutliers --outFileCorMatrix pearson_cor.tab
#plotCorrelation -in all.npz -c pearson -p heatmap -o cor_heatmap.pdf --skipZeros --removeOutliers --colorMap RdYlBu --plotNumbers

#multiBigwigSummary bins -b ../../bigwig/H0_H3k9ac_R?.bigWig -o H0_H3k9ac.npz -bs 200 --smartLabels -p 48
plotCorrelation -in H0_H3k9ac.npz -c pearson -p scatterplot -o H0_H3k9ac_scatter.svg --skipZeros --removeOutliers --plotFileFormat svg

#multiBigwigSummary bins -b ../../bigwig/H24_H3k9ac_R?.bigWig -o H24_H3k9ac.npz -bs 200 --smartLabels -p 48
plotCorrelation -in H24_H3k9ac.npz -c pearson -p scatterplot -o H24_H3k9ac_scatter.svg --skipZeros --removeOutliers --plotFileFormat svg

#multiBigwigSummary bins -b ../../bigwig/H0_H3k27ac_R?.bigWig -o H0_H3k27ac.npz -bs 200 --smartLabels -p 48
plotCorrelation -in H0_H3k27ac.npz -c pearson -p scatterplot -o H0_H3k27ac_scatter.svg --skipZeros --removeOutliers --plotFileFormat svg

#multiBigwigSummary bins -b ../../bigwig/H24_H3k27ac_R?.bigWig -o H24_H3k27ac.npz -bs 200 --smartLabels -p 48
plotCorrelation -in H24_H3k27ac.npz -c pearson -p scatterplot -o H24_H3k27ac_scatter.svg --skipZeros --removeOutliers --plotFileFormat svg

#multiBigwigSummary bins -b ../../bigwig/H0_H3k4m2_R?.bigWig -o H0_H3k4m2.npz -bs 200 --smartLabels -p 48
plotCorrelation -in H0_H3k4m2.npz -c pearson -p scatterplot -o H0_H3k4m2_scatter.svg --skipZeros --removeOutliers --plotFileFormat svg

#multiBigwigSummary bins -b ../../bigwig/H24_H3k4m2_R?.bigWig -o H24_H3k4m2.npz -bs 200 --smartLabels -p 48
plotCorrelation -in H24_H3k4m2.npz -c pearson -p scatterplot -o H24_H3k4m2_scatter.svg --skipZeros --removeOutliers --plotFileFormat svg

#multiBigwigSummary bins -b ../../bigwig/H0_H3kcr_R?.bigWig -o H0_H3kcr.npz -bs 200 --smartLabels -p 48
plotCorrelation -in H0_H3kcr.npz -c pearson -p scatterplot -o H0_H3kcr_scatter.svg --skipZeros --removeOutliers --plotFileFormat svg

#multiBigwigSummary bins -b ../../bigwig/H24_H3kcr_R?.bigWig -o H24_H3kcr.npz -bs 200 --smartLabels -p 48
plotCorrelation -in H24_H3kcr.npz -c pearson -p scatterplot -o H24_H3kcr_scatter.svg --skipZeros --removeOutliers --plotFileFormat svg

