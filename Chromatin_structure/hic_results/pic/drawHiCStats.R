#!/usr/bin/env Rscript

setwd('./')

set.seed(123456)
iter <- 1000

library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

library(reshape2)

data1 <- read.delim("valids.txt", header=T)
df1 <- melt(data1)

data2 <- read.delim("contact_types.txt", header=T)
df2 <- melt(data2)

pdf("NannoEpi_HiC_stats.pdf", useDingbats=FALSE, width=8, height=4, family="Arial")

p1 <- ggbarplot(df1, x="Sample", y="value", fill="variable", palette="npg") +
      scale_x_discrete(limits = rev(levels(df1$Sample))) + coord_flip() +
  theme_bw() + theme(text=element_text(size=14, family = "Arial"), 
        legend.position="top", legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Sample") + ylab("Percent")

p2 <- ggbarplot(df2, x="Sample", y="value", fill="variable", palette="aaas") + 
      scale_x_discrete(limits = rev(levels(df2$Sample))) + coord_flip() +
  theme_bw() + theme(text=element_text(size=14, family = "Arial"),
        legend.position="top", legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color="black")) +
  xlab("") + ylab("Percent")

ggarrange(p1,p2,nrow=1)

dev.off()

