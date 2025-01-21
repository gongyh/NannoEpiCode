#!/usr/bin/env Rscript

setwd('./')

set.seed(123456)

library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

ann1 <- read.delim("H3k9ac/H3k9ac.consensus_peaks.Nearest_PromoterID.txt", header=F, stringsAsFactors=F)
colnames(ann1) <- c("Peak","Gene")

de_peak4 <- read.delim("H3k9ac/deseq2/H0_H3k9acvsH24_H3k9ac/H0_H3k9acvsH24_H3k9ac.deseq2.FDR0.01.results.bed", header=F, stringsAsFactors=F)[, c("V4","V5")]
colnames(de_peak4) <- c("Peak","plogFC")
de_peak4$plogFC <- -1*de_peak4$plogFC

de_gene <- read.delim("VLC24vsH0_logFC.txt", header=T, stringsAsFactors=F)

d4 <- merge(ann1, de_peak4, by="Peak")
d4 <- merge(d4, de_gene, by="Gene")

d4$Mutant <- "H3K9ac"

ann2 <- read.delim("H3k27ac/H3k27ac.consensus_peaks.Nearest_PromoterID.txt", header=F, stringsAsFactors=F)
colnames(ann2) <- c("Peak","Gene")

de_peak6 <- read.delim("H3k27ac/deseq2/H0_H3k27acvsH24_H3k27ac/H0_H3k27acvsH24_H3k27ac.deseq2.FDR0.01.results.bed", header=F, stringsAsFactors=F) [, c("V4","V5")]
colnames(de_peak6) <- c("Peak","plogFC")
de_peak6$plogFC <- -1*de_peak6$plogFC

d6 <- merge(ann2, de_peak6, by="Peak")
d6 <- merge(d6, de_gene, by="Gene")

d6$Mutant <- "H3K27ac"

ann3 <- read.delim("H3kcr/H3kcr.consensus_peaks.Nearest_PromoterID.txt", header=F, stringsAsFactors=F)
colnames(ann3) <- c("Peak","Gene")

de_peak8 <- read.delim("H3kcr/deseq2/H0_H3kcrvsH24_H3kcr/H0_H3kcrvsH24_H3kcr.deseq2.FDR0.01.results.bed", header=F, stringsAsFactors=F) [, c("V4","V5")]
colnames(de_peak8) <- c("Peak","plogFC")
de_peak8$plogFC <- -1*de_peak8$plogFC

d8 <- merge(ann3, de_peak8, by="Peak")
d8 <- merge(d8, de_gene, by="Gene")

d8$Mutant <- "Kcr"

ann4 <- read.delim("H3k4m2/H3k4m2.consensus_peaks.Nearest_PromoterID.txt", header=F, stringsAsFactors=F)
colnames(ann4) <- c("Peak","Gene")

de_peak9 <- read.delim("H3k4m2/deseq2/H0_H3k4m2vsH24_H3k4m2/H0_H3k4m2vsH24_H3k4m2.deseq2.FDR0.01.results.bed", header=F, stringsAsFactors=F) [, c("V4","V5")]
colnames(de_peak9) <- c("Peak","plogFC")
de_peak9$plogFC <- -1*de_peak9$plogFC

d9 <- merge(ann3, de_peak9, by="Peak")
d9 <- merge(d9, de_gene, by="Gene")

d9$Mutant <- "H3K4me2"

data <- rbind(d4,d6,d8,d9)

pdf("NannoEpi_Scatter_oneline.pdf", useDingbats=FALSE, width=11.5, height=3, family="Arial")

data$Mutant <- factor(data$Mutant, ordered=TRUE, levels = c("H3K9ac", "H3K27ac", "Kcr", "H3K4me2"))
str(data)

ggscatter(data, x="plogFC", y="logFC", color="black", shape=21, size=0.3, facet.by="Mutant",nrow=1,
          add="reg.line",  add.params=list(color = "blue", fill = "lightgray"),
          conf.int = TRUE, cor.coef = TRUE, 
          cor.coeff.args = list(method = "pearson", label.x = -3, label.sep = "\n")) +
    theme_bw() + theme(text=element_text(size=14, family = "Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Peak log2FC") + ylab("Gene log2FC")

dev.off()

