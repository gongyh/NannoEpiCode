#!/usr/bin/env Rscript

setwd('./')

df_5mc <- read.delim('genePromoter_5mC_diff.txt',header=F,stringsAsFactors=F)
colnames(df_5mc) <- c("Gene","Diff_5mC")

df_deg <- read.delim('VLC24vsH0_logFC.txt',header=T,stringsAsFactors=F)
df_deg2 <- read.delim('H0_H24_DEG.txt',header=T,stringsAsFactors=F)

df <- merge(df_5mc, merge(df_deg, df_deg2))

library(ggpubr)
library(reshape2)

library(extrafont)
loadfonts(device = "pdf")

pdf("5mC_log2fc.pdf", width=6, height=6, family="Arial", useDingbats = FALSE)

ggscatterhist(df, x="logFC", y="Diff_5mC", color="DEG", palette=c("blue","gray","red"), size=0.5, ylim=c(-1,1),
              margin.params = list(fill = "DEG", color = "black", size = 0.2),
              xlab="log2fc of genes", ylab="differences of 5mC levels (%)") + theme_bw() +
  theme(text=element_text(size=16,  family="Arial"),legend.position="none",
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

comparisons <- list(c("Up","Not"),c("Not","Down"),c("Up","Down"))
ggboxplot(df, x="DEG",y="Diff_5mC",fill="DEG",palette=c("red","gray","blue"), width=0.6, 
          order=c("Up","Not","Down"), outlier.shape = NA) + 
  coord_cartesian(ylim = c(-0.3,0.3)) +
  xlab("DEGs") + ylab("differences of 5mC levels (%)") + theme_bw() +
  theme(text=element_text(size=16,  family="Arial"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), comparisons = comparisons, 
                     method = "t.test", label.y=c(0.22,0.25,0.28), tip.length=0.003)

dev.off()

