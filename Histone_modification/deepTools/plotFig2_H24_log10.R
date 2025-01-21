#!/usr/bin/env Rscript

library(reshape2)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

library(stringr)

options(scipen=1e6)

abd <- read.delim("HC_VLC.TMM.TPM.matrix")
abd$Control <- abd$Control + 1
abd$VLC <- abd$VLC + 1
#str(abd)

H24_H3k9ac_k5 <- read.delim("H24_H3k9ac_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
t2g <- function(data){str_replace(data,"\\..*","")}
H24_H3k9ac_gids <- sapply(H24_H3k9ac_k5$name, t2g)
H24_H3k9ac_group <- data.frame(Genes=H24_H3k9ac_gids,Group=H24_H3k9ac_k5$deepTools_group)
H24_H3k9ac_df <- merge(abd,H24_H3k9ac_group)

H24_H3k27ac_k5 <- read.delim("H24_H3k27ac_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
H24_H3k27ac_gids <- sapply(H24_H3k27ac_k5$name, t2g)
H24_H3k27ac_group <- data.frame(Genes=H24_H3k27ac_gids,Group=H24_H3k27ac_k5$deepTools_group)
H24_H3k27ac_df <- merge(abd,H24_H3k27ac_group)

H24_H3kcr_k5 <- read.delim("H24_H3kcr_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
H24_H3kcr_gids <- sapply(H24_H3kcr_k5$name, t2g)
H24_H3kcr_group <- data.frame(Genes=H24_H3kcr_gids,Group=H24_H3kcr_k5$deepTools_group)
H24_H3kcr_df <- merge(abd,H24_H3kcr_group)

H24_H3k4m2_k5 <- read.delim("H24_H3k4m2_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
H24_H3k4m2_gids <- sapply(H24_H3k4m2_k5$name, t2g)
H24_H3k4m2_group <- data.frame(Genes=H24_H3k4m2_gids,Group=H24_H3k4m2_k5$deepTools_group)
H24_H3k4m2_df <- merge(abd,H24_H3k4m2_group)

pdf("NannoEpi_Fig2_H24_log10.pdf", useDingbats = FALSE, width = 8, height = 6, family="Arial")

comparisons <- list(c("cluster_4","cluster_5"),c("cluster_3","cluster_4"),
                    c("cluster_2","cluster_3"),c("cluster_1","cluster_2"))

p1 <- ggboxplot(H24_H3k9ac_df,x="Group",y="VLC",fill="Group",palette="npg", width=0.5, add="none", bxp.errorbar=T, 
                size=0.2, outlier.size=0.05, yscale="log10", format.scale = TRUE) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Clusters of H3k9ac (H24)") + ylab("Gene expression (TPM+1, log10)") + 
  stat_compare_means(aes(label = ..p.signif..), method="wilcox.test", comparisons=comparisons, family="Arial")

p2 <- ggboxplot(H24_H3k27ac_df,x="Group",y="VLC",fill="Group",palette="npg", width=0.5,
          add="none", bxp.errorbar=T, size=0.2, outlier.size=0.05, yscale="log10", format.scale = TRUE) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Clusters of H3k27ac (H24)") + ylab("Gene expression (TPM+1, log10)") + 
  stat_compare_means(aes(label = ..p.signif..), method="wilcox.test", comparisons=comparisons, family="Arial")

p3 <- ggboxplot(H24_H3kcr_df,x="Group",y="VLC",fill="Group",palette="npg", width=0.5,
          add="none", bxp.errorbar=T, size=0.2, outlier.size=0.05, yscale="log10", format.scale = TRUE) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Clusters of H3kcr (H24)") + ylab("Gene expression (TPM+1, log10)") + 
  stat_compare_means(aes(label = ..p.signif..), method="wilcox.test", comparisons=comparisons, family="Arial")

p4 <- ggboxplot(H24_H3k4m2_df,x="Group",y="VLC",fill="Group",palette="npg", width=0.5,
          add="none", bxp.errorbar=T, size=0.2, outlier.size=0.05, yscale="log10", format.scale = TRUE) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Clusters of H3k4m2 (H24)") + ylab("Gene expression (TPM+1, log10)") + 
  stat_compare_means(aes(label = ..p.signif..), method="wilcox.test", comparisons=comparisons, family="Arial")

ggarrange(p1,p2,p3,p4,ncol=2,nrow=2,legend="none")

dev.off()

