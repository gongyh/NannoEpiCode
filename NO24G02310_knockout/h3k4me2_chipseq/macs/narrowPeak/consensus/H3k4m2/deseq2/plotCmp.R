#!/usr/bin/env Rscript

setwd('./')

set.seed(123456)

library(dplyr)
library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

library(pheatmap)
library(ComplexHeatmap)
library(circlize)

data <- read.delim("H3k4m2.consensus_peaks.sample.dists_mqc.tsv", comment.char="#", header=T, stringsAsFactors=F, row.names=1)

str(data)

pdf("mut_cmp.pdf", useDingbats=FALSE, width=5, height=4, family="Arial")

Heatmap(data, border = "black", cluster_rows = F, cluster_columns = F, column_names_rot = 30)

dev.off()

