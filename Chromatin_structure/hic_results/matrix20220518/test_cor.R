#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

HC_df <- read.table('HC_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)
LC_df <- read.table('LC_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)

colnames(HC_df) <- c("chrom", "start", "end", "E1")
colnames(LC_df) <- c("chrom", "start", "end", "E1")

HC_df$group <- "HC"
LC_df$group <- "LC"

bin_gd <- read.table('bins_40000_gd.txt', sep="\t", header=T)
bin_gc <- read.table('genome_gc_40000.txt', sep="\t", header=T)
bin_at <- bin_gc
bin_at$AT <- 1 - bin_at$GC

data_HC <- merge(HC_df, bin_gd)
data_LC <- merge(LC_df, bin_gd)

data_HC$chr <- as.numeric(sub("chr","",data_HC$chrom))
data_LC$chr <- as.numeric(sub("chr","",data_LC$chrom))

data2_HC <- merge(HC_df, bin_at)
data2_LC <- merge(LC_df, bin_at)

data2_HC$chr <- as.numeric(sub("chr","",data2_HC$chrom))
data2_LC$chr <- as.numeric(sub("chr","",data2_LC$chrom))

pdf("NannoEpi_HiC_AB_E1_GD.pdf", useDingbats=FALSE, width=9, height=10, family="Arial", onefile=F)

p1 <- ggscatter(data_HC, x="GD", y="E1", size=1, add = "reg.line", conf.int = TRUE, facet.by="chr", ncol=8,
      cor.coef = TRUE, cor.coeff.args = list(method = "pearson", color="red", label.sep = "\n")) +
  theme_bw() +
  theme(text=element_text(size=8, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

p2 <- ggscatter(data_LC, x="GD", y="E1", size=1, add = "reg.line", conf.int = TRUE, facet.by="chr", ncol=8,
      cor.coef = TRUE, cor.coeff.args = list(method = "pearson", color="red", label.sep = "\n")) +
  theme_bw() +
  theme(text=element_text(size=8, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

ggarrange(p1,p2,ncol=1)

dev.off()

pdf("NannoEpi_HiC_AB_E1_AT.pdf", useDingbats=FALSE, width=9, height=10, family="Arial", onefile=F)

p3 <- ggscatter(data2_HC, x="AT", y="E1", size=1, add = "reg.line", conf.int = TRUE, facet.by="chr", ncol=8,
      cor.coef = TRUE, cor.coeff.args = list(method = "pearson", color="red", label.sep = "\n")) +
  theme_bw() +
  theme(text=element_text(size=8, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

p4 <- ggscatter(data2_LC, x="AT", y="E1", size=1, add = "reg.line", conf.int = TRUE, facet.by="chr", ncol=8,
      cor.coef = TRUE, cor.coeff.args = list(method = "pearson", color="red", label.sep = "\n")) +
  theme_bw() +
  theme(text=element_text(size=8, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

ggarrange(p3,p4,ncol=1)

dev.off()
