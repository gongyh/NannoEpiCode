#!/usr/bin/env Rscript

setwd('./')

library(GenomicFeatures)

txdb <- makeTxDbFromGFF(file='genes.gtf')

library(ChIPseeker)
library(ggupset)
library(ggimage)

promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)


peak1 <- readPeakFile("H0_H3k4m2_R1_summits.bed")
#covplot(peak1, weightCol="V5")

#tagMatrix <- getTagMatrix(peak1, windows=promoter)
#plotAvgProf(tagMatrix, xlim=c(-2000, 2000))

peak1Anno <- annotatePeak(peak1, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
# plotAnnoPie(peak1Anno)
# plotAnnoBar(peak1Anno)
# upsetplot(peak1Anno, vennpie=TRUE)
# plotDistToTSS(peak1Anno)

p1df <- as.data.frame(peak1Anno)
p1f <- p1df[(p1df$distanceToTSS>0)&(p1df$distanceToTSS<1500),]

write.table(p1f,"H0_H3K4me2_R1_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak2 <- readPeakFile("H0_H3k4m2_R2_summits.bed")
peak2Anno <- annotatePeak(peak2, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p2df <- as.data.frame(peak2Anno)
p2f <- p2df[(p2df$distanceToTSS>0)&(p2df$distanceToTSS<1500),]

write.table(p2f,"H0_H3K4me2_R2_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak3 <- readPeakFile("H24_H3k4m2_R1_summits.bed")
peak3Anno <- annotatePeak(peak3, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p3df <- as.data.frame(peak3Anno)
p3f <- p3df[(p3df$distanceToTSS>0)&(p3df$distanceToTSS<1500),]

write.table(p3f,"H24_H3K4me2_R1_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak4 <- readPeakFile("H24_H3k4m2_R2_summits.bed")
peak4Anno <- annotatePeak(peak4, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p4df <- as.data.frame(peak4Anno)
p4f <- p4df[(p4df$distanceToTSS>0)&(p4df$distanceToTSS<1500),]

write.table(p4f,"H24_H3K4me2_R2_summits_annotations.txt",sep='\t',quote=F,row.names = F)

