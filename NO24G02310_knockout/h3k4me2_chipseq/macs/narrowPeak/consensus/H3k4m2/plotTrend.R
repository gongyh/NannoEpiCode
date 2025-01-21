#!/usr/bin/env Rscript

setwd('./')

set.seed(123456)

library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

ann <- read.delim("H3k4m2.consensus_peaks.Nearest_PromoterID.txt", header=F, stringsAsFactors=F)
colnames(ann) <- c("Peak","Gene")

de_peak4 <- read.delim("M4_logFC.txt", header=F, stringsAsFactors=F)
colnames(de_peak4) <- c("Peak","plogFC")

de_gene4 <- read.delim("M4vsWT_logFC.txt", header=T, stringsAsFactors=F)

d4 <- merge(ann, de_peak4, by="Peak")
d4 <- merge(d4, de_gene4, by="Gene")

d4$Mutant <- "M4"

de_peak6 <- read.delim("M6_logFC.txt", header=F, stringsAsFactors=F)
colnames(de_peak6) <- c("Peak","plogFC")

de_gene6 <- read.delim("M6vsWT_logFC.txt", header=T, stringsAsFactors=F)

d6 <- merge(ann, de_peak6, by="Peak")
d6 <- merge(d6, de_gene6, by="Gene")

d6$Mutant <- "M6"

data <- rbind(d4,d6)

pdf("NannoEpi_Fig4m2_M4M6_Scatter.pdf", useDingbats=FALSE, width=8, height=4, family="Arial")

ggscatter(data, x="plogFC", y="logFC", color="black", shape=21, size=0.3, facet.by="Mutant",
          add="reg.line",  add.params=list(color = "blue", fill = "lightgray"),
          conf.int = TRUE, cor.coef = TRUE, 
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
    theme_bw() + theme(text=element_text(size=14, family = "Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Peak log2FC") + ylab("Gene log2FC")

dev.off()

