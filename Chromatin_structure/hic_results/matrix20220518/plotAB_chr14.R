#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

HC_R1_df <- read.table('NannoH0_compartments.cis.E1.bedgraph', sep="\t", header=F)
HC_R2_df <- read.table('H0rep_compartments.cis.E1.bedgraph', sep="\t", header=F)
LC_R1_df <- read.table('NannoH24_compartments.cis.E1.bedgraph', sep="\t", header=F)
LC_R2_df <- read.table('H24rep_compartments.cis.E1.bedgraph', sep="\t", header=F)

HC_R1_df$pos <- (HC_R1_df$V2 + HC_R1_df$V3)/2
HC_R2_df$pos <- (HC_R2_df$V2 + HC_R2_df$V3)/2
LC_R1_df$pos <- (LC_R1_df$V2 + LC_R1_df$V3)/2
LC_R2_df$pos <- (LC_R2_df$V2 + LC_R2_df$V3)/2

HC_R1_df$sign <- HC_R1_df$V4 < 0.0
HC_R2_df$sign <- HC_R2_df$V4 < 0.0
LC_R1_df$sign <- LC_R1_df$V4 < 0.0
LC_R2_df$sign <- LC_R2_df$V4 < 0.0

pdf("NannoEpi_HiC_AB_chr14.pdf", useDingbats=FALSE, width=5, height=2, family="Arial", onefile=F)

HC_data <- HC_R1_df[HC_R1_df$V1=="chr14",]
LC_data <- LC_R1_df[LC_R1_df$V1=="chr14",]

HC_data$grp <- 'HC'
LC_data$grp <- 'LC'

data <- rbind(HC_data, LC_data)

ggdotchart(data, x="pos", y="V4", sorting="none", add = "segment", add.params=list(size=0.5,color="sign"),
           facet.by="grp", ncol=1, strip.position='right', size=0.5, color="sign", dot.size=0.25, palette=c("red","blue")) +
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black")) +
  xlab("Chromosome 14") + ylab("")

dev.off()

