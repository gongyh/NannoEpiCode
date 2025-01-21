#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

library(gtools)

HC_df <- read.table('H0_summit_count_bin5k.tsv', sep="\t", header=F)
HC_df$Group <- "HC"

LC_df <- read.table('H24_summit_count_bin5k.tsv', sep="\t", header=F)
LC_df$Group <- "LC"

data <- rbind(HC_df, LC_df)
colnames(data) <- c("Bin","Chr","Start","End","Count","HCLC")

data$Group <- paste0(data$Chr,"_",data$HCLC)

dl <- mixedsort(unique(data$Group), decreasing=T)

data$Group <- factor(data$Group, levels=dl, ordered=T)

pdf("NannoEpi_MNase_density.pdf", useDingbats=FALSE, width=5, height=8, family="Arial", onefile=F)

ggbarplot(data, "Group", y=5, color="Count", sort.val = "none", sort.by.groups = F, orientation = "horiz") + 
  scale_y_continuous(expand=expansion(0,0), position="right") + 
  scale_colour_gradient2(low = "darkgreen", mid="yellow", high = "red", midpoint=13) +
  clean_theme() + ylab("") + xlab("") +
  theme(text=element_text(size=11, family="Arial"), legend.position=c(0.75, 0.2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.x.top = element_line(),
        axis.ticks.x.top = element_line(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

dev.off()

