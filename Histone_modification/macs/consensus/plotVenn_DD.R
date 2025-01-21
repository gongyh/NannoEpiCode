#!/usr/bin/env Rscript

library(Vennerable)
#library(VennDiagram)
library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

colors <- get_palette("npg", 4)

pdf("NannoEpi_VennDD.pdf", useDingbats = FALSE, width = 4, height = 4, family="Arial")

H3k9ac <- read.table('H3k9ac_DEPeaks_DEGs.txt',header=F,stringsAsFactors=F)$V1
H3k27ac <- read.table('H3k27ac_DEPeaks_DEGs.txt',header=F,stringsAsFactors=F)$V1
H3kcr <- read.table('H3kcr_DEPeaks_DEGs.txt',header=F,stringsAsFactors=F)$V1
H3k4m2 <- read.table('H3k4m2_DEPeaks_DEGs.txt',header=F,stringsAsFactors=F)$V1

x = list("H3K9ac"=H3k9ac,"H3K27ac"=H3k27ac,"Kcr"=H3kcr,"H3K4me2"=H3k4m2)
p <- Venn(x)
plot(p, doEuler = TRUE)
dev.off()

