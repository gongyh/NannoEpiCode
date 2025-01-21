#!/usr/bin/env Rscript

library(reshape2)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

library(stringr)

abd <- read.delim("HC_VLC.TMM.TPM.matrix")

H0_H3k9ac_k5 <- read.delim("H0_H3k9ac_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
t2g <- function(data){str_replace(data,"\\..*","")}
H0_H3k9ac_gids <- sapply(H0_H3k9ac_k5$name, t2g)
H0_H3k9ac_group <- data.frame(Genes=H0_H3k9ac_gids,Group=H0_H3k9ac_k5$deepTools_group)
H0_H3k9ac_df <- merge(abd,H0_H3k9ac_group)

H0_H3k27ac_k5 <- read.delim("H0_H3k27ac_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
H0_H3k27ac_gids <- sapply(H0_H3k27ac_k5$name, t2g)
H0_H3k27ac_group <- data.frame(Genes=H0_H3k27ac_gids,Group=H0_H3k27ac_k5$deepTools_group)
H0_H3k27ac_df <- merge(abd,H0_H3k27ac_group)

H0_H3kcr_k5 <- read.delim("H0_H3kcr_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
H0_H3kcr_gids <- sapply(H0_H3kcr_k5$name, t2g)
H0_H3kcr_group <- data.frame(Genes=H0_H3kcr_gids,Group=H0_H3kcr_k5$deepTools_group)
H0_H3kcr_df <- merge(abd,H0_H3kcr_group)

H0_H3k4m2_k5 <- read.delim("H0_H3k4m2_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
H0_H3k4m2_gids <- sapply(H0_H3k4m2_k5$name, t2g)
H0_H3k4m2_group <- data.frame(Genes=H0_H3k4m2_gids,Group=H0_H3k4m2_k5$deepTools_group)
H0_H3k4m2_df <- merge(abd,H0_H3k4m2_group)

pdf("NannoEpi_Fig2_H0.pdf", useDingbats = FALSE, width = 7, height = 6, family="Arial")

comparisons <- list(c("cluster_1","cluster_5"),c("cluster_1","cluster_4"),c("cluster_1","cluster_3"),
c("cluster_1","cluster_2"),c("cluster_2","cluster_3"),c("cluster_3","cluster_4"),c("cluster_4","cluster_5"))

p1 <- ggboxplot(H0_H3k9ac_df,x="Group",y="Control",fill="Group",palette="npg", width=0.5,
          add="none", bxp.errorbar=T, size=0.2, outlier.shape = NA, ylim=c(0,300)) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Clusters of H3k9ac (H0)") + ylab("TPM of genes") + 
  stat_compare_means(method="wilcox.test", comparisons=comparisons, family="Arial", 
        tip.length=0.001, size=3, label.y=c(290,270,250,230,210,190,170))

p2 <- ggboxplot(H0_H3k27ac_df,x="Group",y="Control",fill="Group",palette="npg", width=0.5,
          add="none", bxp.errorbar=T, size=0.2, outlier.shape = NA, ylim=c(0,300)) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Clusters of H3k27ac (H0)") + ylab("TPM of genes") + 
  stat_compare_means(method="wilcox.test", comparisons=comparisons, family="Arial", 
        tip.length=0.001, size=3, label.y=c(290,270,250,230,210,190,170))

p3 <- ggboxplot(H0_H3kcr_df,x="Group",y="Control",fill="Group",palette="npg", width=0.5,
          add="none", bxp.errorbar=T, size=0.2, outlier.shape = NA, ylim=c(0,300)) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Clusters of H3kcr (H0)") + ylab("TPM of genes") + 
  stat_compare_means(method="wilcox.test", comparisons=comparisons, family="Arial", 
        tip.length=0.001, size=3, label.y=c(290,270,250,230,210,190,170))

p4 <- ggboxplot(H0_H3k4m2_df,x="Group",y="Control",fill="Group",palette="npg", width=0.5,
          add="none", bxp.errorbar=T, size=0.2, outlier.shape = NA, ylim=c(0,300)) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Clusters of H3k4m2 (H0)") + ylab("TPM of genes") + 
  stat_compare_means(method="wilcox.test", comparisons=comparisons, family="Arial", 
        tip.length=0.001, size=3, label.y=c(290,270,250,230,210,190,170))

ggarrange(p1,p2,p3,p4,ncol=2,nrow=2,legend="none")

dev.off()

