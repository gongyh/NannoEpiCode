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

data <- read.delim("emissions_5.txt", header=T, stringsAsFactors=F, row.names=1)

pdf("emissions.pdf", useDingbats=FALSE, width=7, height=2, family="Arial")

#pheatmap(data, display_numbers = T, cluster_rows = F, cluster_cols = F, 
#         border_color = "black", legend = F, number_color = "black",
#         color = colorRampPalette(c('white','red'))(100))

col_fun = colorRamp2(c(0, 1), c("white", "red"))

anno <- c("Deficient in all marks", "Enriched in nucleosomes", "Enriched in H3k4m2",
          "Enriched in all the four histones", "Enriched in H3k9ac & H3k27ac & H3kcr")

Heatmap(data, col = col_fun, border = "black", cluster_rows = F, cluster_columns = F, 
        row_names_side = "left", column_names_side = "top", column_names_rot = 0, column_names_centered = T,
        show_heatmap_legend = F, right_annotation = rowAnnotation(ann = anno_text(anno)),
    cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
        grid.text(sprintf("%.2f",data[i, j]), x, y)
    })

dev.off()

