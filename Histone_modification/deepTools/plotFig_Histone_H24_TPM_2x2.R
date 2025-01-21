#!/usr/bin/env Rscript

library(reshape2)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

library(stringr)

options(scipen=1e6)

abd <- read.delim("HC_VLC.TMM.TPM.matrix")
abd$Control <- abd$Control
abd$VLC <- abd$VLC
#str(abd)

H24_H3k9ac_k5 <- read.delim("H24_H3k9ac_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
t2g <- function(data){str_replace(data,"\\..*","")}
H24_H3k9ac_gids <- sapply(H24_H3k9ac_k5$name, t2g)
H24_H3k9ac_group <- data.frame(Genes=H24_H3k9ac_gids,Group=H24_H3k9ac_k5$deepTools_group)
H24_H3k9ac_df <- merge(abd,H24_H3k9ac_group)
H24_H3k9ac_df$Type <- "H3K9ac"

H24_H3k27ac_k5 <- read.delim("H24_H3k27ac_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
H24_H3k27ac_gids <- sapply(H24_H3k27ac_k5$name, t2g)
H24_H3k27ac_group <- data.frame(Genes=H24_H3k27ac_gids,Group=H24_H3k27ac_k5$deepTools_group)
H24_H3k27ac_df <- merge(abd,H24_H3k27ac_group)
H24_H3k27ac_df$Type <- "H3K27ac"

H24_H3kcr_k5 <- read.delim("H24_H3kcr_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
H24_H3kcr_gids <- sapply(H24_H3kcr_k5$name, t2g)
H24_H3kcr_group <- data.frame(Genes=H24_H3kcr_gids,Group=H24_H3kcr_k5$deepTools_group)
H24_H3kcr_df <- merge(abd,H24_H3kcr_group)
H24_H3kcr_df$Type <- "Kcr"

H24_H3k4m2_k5 <- read.delim("H24_H3k4m2_k5.plotHeatmap.sorted.bed", stringsAsFactors=F)
H24_H3k4m2_gids <- sapply(H24_H3k4m2_k5$name, t2g)
H24_H3k4m2_group <- data.frame(Genes=H24_H3k4m2_gids,Group=H24_H3k4m2_k5$deepTools_group)
H24_H3k4m2_df <- merge(abd,H24_H3k4m2_group)
H24_H3k4m2_df$Type <- "H3K4me2"

data <- rbind(H24_H3k9ac_df,H24_H3k27ac_df,H24_H3kcr_df,H24_H3k4m2_df)
data$Type <- factor(data$Type, ordered=TRUE, levels = c("H3K9ac", "H3K27ac", "Kcr", "H3K4me2"))
str(data)

pdf("NannoEpi_Histone_H24_TPM_2x2.pdf", useDingbats = FALSE, width = 5, height = 7, family="Arial")

comparisons <- list(c("cluster_4","cluster_5"),c("cluster_3","cluster_4"),
                    c("cluster_2","cluster_3"),c("cluster_1","cluster_2"))

ggboxplot(data,x="Group",y="Control",fill="Group",palette="npg", width=0.5, add="none", bxp.errorbar=T, 
                size=0.2, outlier.shape=NA, facet.by="Type", nrow=2) + theme_bw() +
  scale_x_discrete(labels = c("cluster_1"="C1","cluster_2"="C2","cluster_3"="C3","cluster_4"="C4","cluster_5"="C5")) +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) + coord_cartesian(ylim=c(0,250)) +
  xlab("") + ylab("Gene expression (TPM)") + 
  stat_compare_means(aes(label = ..p.signif..), method="wilcox.test", comparisons=comparisons, 
                     family="Arial", label.y=c(210,220,230,240), tip.length =0.0005)

dev.off()

