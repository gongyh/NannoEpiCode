#!/usr/bin/Rscript

require(ggplot2)
require(ggpubr)
require(webr)

library(reshape2)
library(extrafont)
#font_import()
loadfonts(device = "pdf")


data <- read.table("macs_annotatePeaks.summary.R1R2.txt", header=T, sep="\t",row.names=1)

str(data)

df <- data.frame(t(data))
df$type <- rownames(df)
dm <- melt(df)

pdf("NannoEpi_anno_pie.pdf", useDingbats = FALSE, width = 4, height = 4, family="Arial")

p1 <- PieDonut(dm[dm$variable=="HC_H3k9ac",], aes(type,count=value), ratioByGroup=F, family="Arial", r0=0.6, r1=0.9, r2=1.2, labelposition=0, pieLabelSize=2) +
      PieDonut(dm[dm$variable=="LC_H3k9ac",], aes(type,count=value), ratioByGroup=F, family="Arial", r0=0.9, r1=1.2, r2=1.2, labelposition=0, pieLabelSize=2)$layers +
      theme(panel.border = element_blank())

p2 <- PieDonut(dm[dm$variable=="HC_H3k27ac",], aes(type,count=value), ratioByGroup=F, family="Arial", r0=0.6, r1=0.9, r2=1.2, labelposition=0, pieLabelSize=2) +
      PieDonut(dm[dm$variable=="LC_H3k27ac",], aes(type,count=value), ratioByGroup=F, family="Arial", r0=0.9, r1=1.2, r2=1.2, labelposition=0, pieLabelSize=2)$layers +
      theme(panel.border = element_blank())

p3 <- PieDonut(dm[dm$variable=="HC_kcr",], aes(type,count=value), ratioByGroup=F, family="Arial", r0=0.6, r1=0.9, r2=1.2, labelposition=0, pieLabelSize=2) +
      PieDonut(dm[dm$variable=="LC_kcr",], aes(type,count=value), ratioByGroup=F, family="Arial", r0=0.9, r1=1.2, r2=1.2, labelposition=0, pieLabelSize=2)$layers +
      theme(panel.border = element_blank())

p4 <- PieDonut(dm[dm$variable=="HC_H3k4m2",], aes(type,count=value), ratioByGroup=F, family="Arial", r0=0.6, r1=0.9, r2=1.2, labelposition=0, pieLabelSize=2) +
      PieDonut(dm[dm$variable=="LC_H3k4m2",], aes(type,count=value), ratioByGroup=F, family="Arial", r0=0.9, r1=1.2, r2=1.2, labelposition=0, pieLabelSize=2)$layers +
      theme(panel.border = element_blank())

ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)

dev.off()

