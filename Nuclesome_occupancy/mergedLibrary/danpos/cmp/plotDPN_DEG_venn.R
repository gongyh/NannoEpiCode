#!/usr/bin/env Rscript

library(VennDiagram)
library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

colors <- get_palette("npg", 4)

pdf("NannoEpi_MNase_venn.pdf", useDingbats = FALSE, width = 8, height = 8, family="Arial")

draw.triple.venn(area1=1873, area2=1892, area3=1512, n12=332,n23=332,n13=0,n123=0,euler.d=T,
       category = c("DEG_Up", "Genes with DPNs", "DEG_Down"), fill = c(colors[2],colors[1],colors[3]), lwd=1,
       cat.cex = 1.45, cat.dist = 0.05, scaled = TRUE, alpha=0.6, cat.col = c(colors[2],colors[1],colors[3]),
       cat.fontfamily = "Arial", cex = 1.45, fontfamily = "Arial", margin=0.2)

dev.off()

