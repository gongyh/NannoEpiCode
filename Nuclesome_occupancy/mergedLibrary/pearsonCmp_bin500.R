#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

HC_R1_df <- read.table('H0_R1.mLb.clN.tab', sep="\t", header=F)
HC_R2_df <- read.table('H0_R2.mLb.clN.tab', sep="\t", header=F)
LC_R1_df <- read.table('H24_R1.mLb.clN.tab', sep="\t", header=F)
LC_R2_df <- read.table('H24_R2.mLb.clN.tab', sep="\t", header=F)

HC_df <- inner_join(HC_R1_df[,c("V1","V4")],HC_R2_df[,c("V1","V4")],by="V1")
LC_df <- inner_join(LC_R1_df[,c("V1","V4")],LC_R2_df[,c("V1","V4")],by="V1")

HC_data <- data.frame(HC_R1=HC_df$V4.x,HC_R2=HC_df$V4.y)
LC_data <- data.frame(LC_R1=LC_df$V4.x,LC_R2=LC_df$V4.y)

pdf("NannoEpi_MNase_replicate_pearson.pdf", useDingbats=FALSE, width=8.5, height=4, family="Arial", onefile=F)

hc <- ggplot(HC_data, aes(HC_R1, HC_R2), xscale="log10", yscale="log10") + 
  geom_hex(bins=60) + stat_cor(label.sep="\n") + theme_bw() +
  scale_fill_gradient(low="#00AFBB", high="#FC4E07", trans="log10") +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("HC rep1") + ylab("HC rep2")

lc <- ggplot(LC_data, aes(LC_R1, LC_R2), xscale="log10", yscale="log10") + 
  geom_hex(bins=60) + stat_cor(label.sep="\n") + theme_bw() +
  scale_fill_gradient(low="#00AFBB", high="#FC4E07", trans="log10") +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("LC rep1") + ylab("LC rep2")

ggarrange(hc, lc, nrow=1, legend="right", common.legend=T)

dev.off()

