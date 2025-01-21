#!/usr/bin/env Rscript

library(VennDiagram)
library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

colors <- get_palette("npg", 4)

pdf("NannoEpi_Fig1C.pdf", useDingbats = FALSE, width = 4, height = 4, family="Arial")

p1 <- draw.pairwise.venn(area1=6228, area2=6130, cross.area = 6055,
       category = c("H3k9ac", "H3k27ac"), fill = c(colors[1],colors[2]), lwd=1,
       cat.cex = 1.45, cat.dist = 0.05, scaled = TRUE, alpha=0.6, cat.col = c(colors[1],colors[2]),
       cat.fontfamily = "Arial", cex = 1.45, fontfamily = "Arial", margin=0.2)

grid.draw(p1)
grid.newpage()

p2 <- draw.pairwise.venn(area1=6228, area2=6202, cross.area = 6125,
       category = c("H3k9ac", "H3kcr"), fill = c(colors[1],colors[3]), lwd=1,
       cat.cex = 1.45, cat.dist = 0.05, scaled = TRUE, alpha=0.6, cat.col = c(colors[1],colors[3]),
       cat.fontfamily = "Arial", cex = 1.45, fontfamily = "Arial", margin=0.2)

grid.draw(p2)
grid.newpage()

p3 <- draw.pairwise.venn(area1=6130, area2=6202, cross.area = 6107,
       category = c("H3k27ac", "H3kcr"), fill = c(colors[2],colors[3]), lwd=1,
       cat.cex = 1.45, cat.dist = 0.05, scaled = TRUE, alpha=0.6, cat.col = c(colors[2],colors[3]),
       cat.fontfamily = "Arial", cex = 1.45, fontfamily = "Arial", margin=0.2)

grid.draw(p3)

H3k9ac <- read.table('H3k9ac_genes.txt',header=F,stringsAsFactors=F)$V1
H3k27ac <- read.table('H3k27ac_genes.txt',header=F,stringsAsFactors=F)$V1
H3kcr <- read.table('H3kcr_genes.txt',header=F,stringsAsFactors=F)$V1
H3k4m2 <- read.table('H3k4m2_genes.txt',header=F,stringsAsFactors=F)$V1

x = list("H3k9ac"=H3k9ac,"H3k27ac"=H3k27ac,"H3kcr"=H3kcr,"H3k4m2"=H3k4m2)
#overlap <- calculate.overlap(x)
grid.newpage()
p4 <- venn.diagram(x, filename=NULL, lwd=1, fill=colors, alpha=0.6, cex = 1, fontfamily = "Arial", cat.cex = 1, 
                   cat.dist = 0.1, cat.col = colors, cat.fontfamily = "Arial", margin=0.2)
grid.draw(p4)

dev.off()

