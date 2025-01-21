#!/usr/bin/env Rscript

library(reshape2)
library(ggpubr)
library(ggsignif)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

library(stringr)
library(dplyr)

options(scipen=1e6)

abd <- read.delim("HC_VLC.TMM.TPM.matrix.noCpMt", stringsAsFactors=F, header=T)
abd$Control <- abd$Control
abd$VLC <- abd$VLC

H0_H3k9ac_genes <- read.delim("H0_H3k9ac_peaks_genes.txt", stringsAsFactors=F, header=F)
H24_H3k9ac_genes <- read.delim("H24_H3k9ac_peaks_genes.txt", stringsAsFactors=F, header=F)
H0_H3k27ac_genes <- read.delim("H0_H3k27ac_peaks_genes.txt", stringsAsFactors=F, header=F)
H24_H3k27ac_genes <- read.delim("H24_H3k27ac_peaks_genes.txt", stringsAsFactors=F, header=F)
H0_H3kcr_genes <- read.delim("H0_H3kcr_peaks_genes.txt", stringsAsFactors=F, header=F)
H24_H3kcr_genes <- read.delim("H24_H3kcr_peaks_genes.txt", stringsAsFactors=F, header=F)
H0_H3k4m2_genes <- read.delim("H0_H3k4m2_peaks_genes.txt", stringsAsFactors=F, header=F)
H24_H3k4m2_genes <- read.delim("H24_H3k4m2_peaks_genes.txt", stringsAsFactors=F, header=F)

H0_H3k9ac_peaks_tpm <- filter(abd, Genes %in% intersect(abd$Genes,H0_H3k9ac_genes$V1))$Control
H0_H3k9ac_nopeaks_tpm <- filter(abd, Genes %in% setdiff(abd$Genes,H0_H3k9ac_genes$V1))$Control
data1 <- data.frame(histone="HC_H3k9ac",TPM=H0_H3k9ac_peaks_tpm,Peaks="Yes")
data2 <- data.frame(histone="HC_H3k9ac",TPM=H0_H3k9ac_nopeaks_tpm,Peaks="No")
wilcox.test(H0_H3k9ac_peaks_tpm,H0_H3k9ac_nopeaks_tpm)
H24_H3k9ac_peaks_tpm <- filter(abd, Genes %in% intersect(abd$Genes,H24_H3k9ac_genes$V1))$VLC
H24_H3k9ac_nopeaks_tpm <- filter(abd, Genes %in% setdiff(abd$Genes,H24_H3k9ac_genes$V1))$VLC
data3 <- data.frame(histone="LC_H3k9ac",TPM=H24_H3k9ac_peaks_tpm,Peaks="Yes")
data4 <- data.frame(histone="LC_H3k9ac",TPM=H24_H3k9ac_nopeaks_tpm,Peaks="No")
wilcox.test(H24_H3k9ac_peaks_tpm,H24_H3k9ac_nopeaks_tpm)

H0_H3k27ac_peaks_tpm <- filter(abd, Genes %in% intersect(abd$Genes,H0_H3k27ac_genes$V1))$Control
H0_H3k27ac_nopeaks_tpm <- filter(abd, Genes %in% setdiff(abd$Genes,H0_H3k27ac_genes$V1))$Control
data5 <- data.frame(histone="HC_H3k27ac",TPM=H0_H3k27ac_peaks_tpm,Peaks="Yes")
data6 <- data.frame(histone="HC_H3k27ac",TPM=H0_H3k27ac_nopeaks_tpm,Peaks="No")
wilcox.test(H0_H3k27ac_peaks_tpm,H0_H3k27ac_nopeaks_tpm)
H24_H3k27ac_peaks_tpm <- filter(abd, Genes %in% intersect(abd$Genes,H24_H3k27ac_genes$V1))$VLC
H24_H3k27ac_nopeaks_tpm <- filter(abd, Genes %in% setdiff(abd$Genes,H24_H3k27ac_genes$V1))$VLC
data7 <- data.frame(histone="LC_H3k27ac",TPM=H24_H3k27ac_peaks_tpm,Peaks="Yes")
data8 <- data.frame(histone="LC_H3k27ac",TPM=H24_H3k27ac_nopeaks_tpm,Peaks="No")
wilcox.test(H24_H3k27ac_peaks_tpm,H24_H3k27ac_nopeaks_tpm)

H0_H3kcr_peaks_tpm <- filter(abd, Genes %in% intersect(abd$Genes,H0_H3kcr_genes$V1))$Control
H0_H3kcr_nopeaks_tpm <- filter(abd, Genes %in% setdiff(abd$Genes,H0_H3kcr_genes$V1))$Control
data9 <- data.frame(histone="HC_H3kcr",TPM=H0_H3kcr_peaks_tpm,Peaks="Yes")
data10 <- data.frame(histone="HC_H3kcr",TPM=H0_H3kcr_nopeaks_tpm,Peaks="No")
wilcox.test(H0_H3kcr_peaks_tpm,H0_H3kcr_nopeaks_tpm)
H24_H3kcr_peaks_tpm <- filter(abd, Genes %in% intersect(abd$Genes,H24_H3kcr_genes$V1))$VLC
H24_H3kcr_nopeaks_tpm <- filter(abd, Genes %in% setdiff(abd$Genes,H24_H3kcr_genes$V1))$VLC
data11 <- data.frame(histone="LC_H3kcr",TPM=H24_H3kcr_peaks_tpm,Peaks="Yes")
data12 <- data.frame(histone="LC_H3kcr",TPM=H24_H3kcr_nopeaks_tpm,Peaks="No")
wilcox.test(H24_H3kcr_peaks_tpm,H24_H3kcr_nopeaks_tpm)

H0_H3k4m2_peaks_tpm <- filter(abd, Genes %in% intersect(abd$Genes,H0_H3k4m2_genes$V1))$Control
H0_H3k4m2_nopeaks_tpm <- filter(abd, Genes %in% setdiff(abd$Genes,H0_H3k4m2_genes$V1))$Control
data13 <- data.frame(histone="HC_H3k4m2",TPM=H0_H3k4m2_peaks_tpm,Peaks="Yes")
data14 <- data.frame(histone="HC_H3k4m2",TPM=H0_H3k4m2_nopeaks_tpm,Peaks="No")
wilcox.test(H0_H3k4m2_peaks_tpm,H0_H3k4m2_nopeaks_tpm)
H24_H3k4m2_peaks_tpm <- filter(abd, Genes %in% intersect(abd$Genes,H24_H3k4m2_genes$V1))$VLC
H24_H3k4m2_nopeaks_tpm <- filter(abd, Genes %in% setdiff(abd$Genes,H24_H3k4m2_genes$V1))$VLC
data15 <- data.frame(histone="LC_H3k4m2",TPM=H24_H3k4m2_peaks_tpm,Peaks="Yes")
data16 <- data.frame(histone="LC_H3k4m2",TPM=H24_H3k4m2_nopeaks_tpm,Peaks="No")
wilcox.test(H24_H3k4m2_peaks_tpm,H24_H3k4m2_nopeaks_tpm)

data <- rbind(data1,data2,data3,data4,data5,data6,data7,data8,data9,
              data10,data11,data12,data13,data14,data15,data16)

stat.test <- compare_means(TPM ~ Peaks, data=data, group.by="histone", method = "wilcox.test")

pdf("NannoEpi_Fig3_tpm.pdf", useDingbats = FALSE, width = 8, height = 4, family="Arial")

ggboxplot(data,x="histone",y="TPM",fill="Peaks",palette="npg", width=0.5, add="none", bxp.errorbar=T,
                size=0.2, outlier.shape=NA) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) + coord_cartesian(ylim=c(0,220)) +
  xlab("") + ylab("TPM") +
  stat_pvalue_manual(
    stat.test, x = "histone", y.position = 200,
    label = "p.signif", position = position_dodge(0.8)
  )
#  geom_signif(data = pdata, aes(xmin=start,y_position=y,xmax=end,annotations=label), manual = TRUE)

dev.off()

