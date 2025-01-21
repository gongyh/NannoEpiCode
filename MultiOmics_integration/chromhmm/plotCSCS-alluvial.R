#!/usr/bin/env Rscript

setwd('./')

set.seed(123456)

library(dplyr)
library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

library(reshape2)
library(ggalluvial)

df <- read.delim("CS_CS.txt", header=T, stringsAsFactors=F, row.names=1)
df$HC <- rownames(df)
data <- melt(df, variable.name="LC", value.name="Genes")
data$Change <- data$HC!=data$LC

pdf("CSCS-alluvial.pdf", useDingbats=FALSE, width=2, height=7, family="Arial")

ggplot(data = data, aes(axis1 = HC, axis2 = LC, y = Genes, fill=HC)) +
  scale_x_discrete(limits = c("HC", "LC"), expand = c(.2, .05)) + xlab("") +
  geom_alluvium(aes(fill = HC)) +
  geom_stratum() + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() + theme(legend.position = "none", axis.text.x = element_text(color="black"))

dev.off()

