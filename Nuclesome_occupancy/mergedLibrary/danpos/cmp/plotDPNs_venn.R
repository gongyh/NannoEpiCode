#!/usr/bin/env Rscript

library(VennDiagram)
library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

colors <- get_palette("npg", 9)

pdf("NannoEpi_MNase_DPNs.pdf", useDingbats = FALSE, width = 8, height = 8, family="Arial")

draw.triple.venn(area1=2031, area2=2030, area3=3201, n12=2005,n23=755,n13=753,n123=747,euler.d=T,
       category = c("Summit", "Point", "Fuzziness"), fill = c(colors[4],colors[5],colors[9]), lwd=1,
       cat.cex = 1.45, cat.dist = 0.05, scaled = TRUE, alpha=0.6, cat.col = c(colors[4],colors[5],colors[9]),
       cat.fontfamily = "Arial", cex = 1.45, fontfamily = "Arial", margin=0.2)

dev.off()

