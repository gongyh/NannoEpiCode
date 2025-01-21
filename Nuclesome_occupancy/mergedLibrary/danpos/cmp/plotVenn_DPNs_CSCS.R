#!/usr/bin/env Rscript

library(VennDiagram)
library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

colors <- get_palette("npg", 4)

pdf("NannoEpi_Venn_DPNs_CSCS.pdf", useDingbats = FALSE, width = 4, height = 4, family="Arial")

p1 <- draw.pairwise.venn(area1=24, area2=29, cross.area = 13,
       category = c("CS5->CS4", "DPNs_down_DEGs"), fill = c(colors[1],colors[2]), lwd=1,
       cat.cex = 1.45, cat.dist = 0.05, scaled = TRUE, alpha=0.6, cat.col = c(colors[1],colors[2]),
       cat.fontfamily = "Arial", cex = 1.45, fontfamily = "Arial", margin=0.2)

dev.off()

