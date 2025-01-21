#!/bin/bash

## plot fingerprint

plotFingerprint --bamfiles ../../H*_H3k27ac_R{1,2}.mLb.clN.sorted.bam ../../Input_H*_R{1,2}.mLb.clN.sorted.bam \
        --plotFile H3k27ac.plotFingerprint.pdf \
        --labels H0_H3k27ac_R1 H0_H3k27ac_R2 H24_H3k27ac_R1 H24_H3k27ac_R2 Input_H0_R1 Input_H0_R2 Input_H24_R1 Input_H24_R2 \
        --outRawCounts H3k27ac.plotFingerprint.raw.txt \
        --outQualityMetrics H3k27ac.plotFingerprint.qcmetrics.txt \
        --skipZeros --numberOfProcessors 48

plotFingerprint --bamfiles ../../H*_H3k9ac_R{1,2}.mLb.clN.sorted.bam ../../Input_H*_R{1,2}.mLb.clN.sorted.bam \
        --plotFile H3k9ac.plotFingerprint.pdf \
        --labels H0_H3k9ac_R1 H0_H3k9ac_R2 H24_H3k9ac_R1 H24_H3k9ac_R2 Input_H0_R1 Input_H0_R2 Input_H24_R1 Input_H24_R2 \
        --outRawCounts H3k9ac.plotFingerprint.raw.txt \
        --outQualityMetrics H3k9ac.plotFingerprint.qcmetrics.txt \
        --skipZeros --numberOfProcessors 48

plotFingerprint --bamfiles ../../H*_H3k4m2_R{1,2}.mLb.clN.sorted.bam ../../Input_H*_R{1,2}.mLb.clN.sorted.bam \
        --plotFile H3k4m2.plotFingerprint.pdf \
        --labels H0_H3k4m2_R1 H0_H3k4m2_R2 H24_H3k4m2_R1 H24_H3k4m2_R2 Input_H0_R1 Input_H0_R2 Input_H24_R1 Input_H24_R2 \
        --outRawCounts H3k4m2.plotFingerprint.raw.txt \
        --outQualityMetrics H3k4m2.plotFingerprint.qcmetrics.txt \
        --skipZeros --numberOfProcessors 48

plotFingerprint --bamfiles ../../H*_H3kcr_R{1,2}.mLb.clN.sorted.bam ../../Input_H*_R{1,2}.mLb.clN.sorted.bam \
        --plotFile H3kcr.plotFingerprint.pdf \
        --labels H0_H3kcr_R1 H0_H3kcr_R2 H24_H3kcr_R1 H24_H3kcr_R2 Input_H0_R1 Input_H0_R2 Input_H24_R1 Input_H24_R2 \
        --outRawCounts H3kcr.plotFingerprint.raw.txt \
        --outQualityMetrics H3kcr.plotFingerprint.qcmetrics.txt \
        --skipZeros --numberOfProcessors 48

