#!/usr/bin/env Rscript

library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

hc1 <- read.delim("HC_R1_5mC_10k.txt", stringsAsFactors=T, header=F)
hc1$Group <- "HC_R1"
hc2 <- read.delim("HC_R2_5mC_10k.txt", stringsAsFactors=T, header=F)
hc2$Group <- "HC_R2"
hc3 <- read.delim("HC_R3_5mC_10k.txt", stringsAsFactors=T, header=F)
hc3$Group <- "HC_R3"
lc1 <- read.delim("LC_R1_5mC_10k.txt", stringsAsFactors=T, header=F)
lc1$Group <- "LC_R1"
lc2 <- read.delim("LC_R2_5mC_10k.txt", stringsAsFactors=T, header=F)
lc2$Group <- "LC_R2"
lc3 <- read.delim("LC_R3_5mC_10k.txt", stringsAsFactors=T, header=F)
lc3$Group <- "LC_R3"

data <- rbind(hc1,hc2,hc3,lc1,lc2,lc3)

pdf("NannoEpi_fig1_5mC.pdf", useDingbats = FALSE, width = 7, height = 8, family="Arial")

ggline(data,x="V1",y="V2",color="Group",palette=c("darkblue","darkblue","darkblue","red","red","red"),
       facet.by="Group",ncol=1,plot_type ="l") + 
  theme_bw() + theme(text=element_text(size=11, family="Arial"), 
        legend.position="none", strip.text = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  scale_x_discrete(breaks = NULL) +
  xlab("Chromosome 30 Mb (resolution: 10Kb)") + ylab("5mC per 10Kb")

dev.off()
