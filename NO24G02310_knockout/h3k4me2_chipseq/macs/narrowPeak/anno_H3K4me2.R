#!/usr/bin/env Rscript

setwd('./')

library(GenomicFeatures)

txdb <- makeTxDbFromGFF(file='genes.gtf')

library(ChIPseeker)
library(ggupset)
library(ggimage)

promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)


peak01 <- readPeakFile("WT_R1_summits.bed")
peak01Anno <- annotatePeak(peak01, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p01df <- as.data.frame(peak01Anno)
p01f <- p01df[(p01df$distanceToTSS>0)&(p01df$distanceToTSS<1500),]
write.table(p01f,"WT_R1_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak02 <- readPeakFile("WT_R2_summits.bed")
peak02Anno <- annotatePeak(peak02, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p02df <- as.data.frame(peak02Anno)
p02f <- p02df[(p02df$distanceToTSS>0)&(p02df$distanceToTSS<1500),]
write.table(p02f,"WT_R2_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak03 <- readPeakFile("WT_R3_summits.bed")
peak03Anno <- annotatePeak(peak03, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p03df <- as.data.frame(peak03Anno)
p03f <- p03df[(p03df$distanceToTSS>0)&(p03df$distanceToTSS<1500),]
write.table(p03f,"WT_R3_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak1 <- readPeakFile("M4_R1_summits.bed")
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

write.table(p1f,"M4_R1_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak2 <- readPeakFile("M4_R2_summits.bed")
peak2Anno <- annotatePeak(peak2, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p2df <- as.data.frame(peak2Anno)
p2f <- p2df[(p2df$distanceToTSS>0)&(p2df$distanceToTSS<1500),]

write.table(p2f,"M4_R2_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak23 <- readPeakFile("M4_R3_summits.bed")
peak23Anno <- annotatePeak(peak23, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p23df <- as.data.frame(peak23Anno)
p23f <- p23df[(p23df$distanceToTSS>0)&(p23df$distanceToTSS<1500),]

write.table(p23f,"M4_R3_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak3 <- readPeakFile("M6_R1_summits.bed")
peak3Anno <- annotatePeak(peak3, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p3df <- as.data.frame(peak3Anno)
p3f <- p3df[(p3df$distanceToTSS>0)&(p3df$distanceToTSS<1500),]

write.table(p3f,"M6_R1_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak4 <- readPeakFile("M6_R2_summits.bed")
peak4Anno <- annotatePeak(peak4, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p4df <- as.data.frame(peak4Anno)
p4f <- p4df[(p4df$distanceToTSS>0)&(p4df$distanceToTSS<1500),]

write.table(p4f,"M6_R2_summits_annotations.txt",sep='\t',quote=F,row.names = F)

peak43 <- readPeakFile("M6_R3_summits.bed")
peak43Anno <- annotatePeak(peak43, tssRegion=c(0, 1500), TxDb=txdb, level="gene", overlap="all")
p43df <- as.data.frame(peak43Anno)
p43f <- p4df[(p43df$distanceToTSS>0)&(p43df$distanceToTSS<1500),]

write.table(p43f,"M6_R3_summits_annotations.txt",sep='\t',quote=F,row.names = F)

