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

df <- read.delim("CS_CS.txt", header=T, stringsAsFactors=F, row.names=1)

data <- as.matrix(df)

pdf("CSCS2.pdf", useDingbats=FALSE, width=3.5, height=2.3, family="Arial")

#pheatmap(data, display_numbers = T, cluster_rows = F, cluster_cols = F, 
#         border_color = "black", legend = F, number_color = "black",
#         color = colorRampPalette(c('white','red'))(100))

data2 <- data/diag(data)
col_fun <- colorRamp2(c(0, 1), c("white", "red"))

Heatmap(data2, col=col_fun, border = "black", cluster_rows = F, cluster_columns = F, 
        row_names_side = "left", column_names_side = "top", column_names_rot = 0, column_names_centered = T,
        show_heatmap_legend = F, row_title = "HC", column_title = "LC",
    cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
        grid.text(sprintf("%d",data[i, j]), x, y)
    })

dev.off()

