#!/usr/bin/env Rscript

setwd('./')

set.seed(123456)
iter <- 1000

library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

data4 <- read.delim("M4_DEGs.txt", header=T, stringsAsFactors=F)
data6 <- read.delim("M6_DEGs.txt", header=T, stringsAsFactors=F)

##################################
H3K4m2_peaks <- read.delim("../M4vsWT_H3K4me2_shift.txt", header=F, stringsAsFactors=F)
num <- 456
H3K4m2_up_ratio <- 11/6

d13 <- data.frame(group="Peaks_DE",ratio=H3K4m2_up_ratio,updn="Far",histone="H3K4me2_M4")

udr <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$V1, num)
  sg_all <- data4[data4$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(ratio, udr)
}

p7 <- sum(udr>H3K4m2_up_ratio)/length(udr)

d14 <- data.frame(group="Random",ratio=udr,updn="Far",histone="H3K4me2_M4")

num <- 6577
H3K4m2_dn_ratio <- 27/174

d15 <- data.frame(group="Peaks_DE",ratio=H3K4m2_dn_ratio,updn="Near",histone="H3K4me2_M4")

udr2 <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$V1, num)
  sg_all <- data4[data4$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg_all[sg_all$DEG=="Up","Gene"])
  udr2 <- c(ratio, udr2)
}

p8 <- sum(udr2>H3K4m2_dn_ratio)/length(udr2)

d16 <- data.frame(group="Random",ratio=udr2,updn="Near",histone="H3K4me2_M4")

##################################
H3K4m2_peaks <- read.delim("../M6vsWT_H3K4me2_shift.txt", header=F, stringsAsFactors=F)
num <- 399
H3K4m2_up_ratio <- 4/8

d17 <- data.frame(group="Peaks_DE",ratio=H3K4m2_up_ratio,updn="Far",histone="H3K4me2_M6")

udr <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$V1, num)
  sg_all <- data6[data6$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(ratio, udr)
}

p7 <- sum(udr>H3K4m2_up_ratio)/length(udr)

d18 <- data.frame(group="Random",ratio=udr,updn="Far",histone="H3K4me2_M6")

num <- 5624
H3K4m2_dn_ratio <- 58/137

d19 <- data.frame(group="Peaks_DE",ratio=H3K4m2_dn_ratio,updn="Near",histone="H3K4me2_M6")

udr2 <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$V1, num)
  sg_all <- data6[data6$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg_all[sg_all$DEG=="Up","Gene"])
  udr2 <- c(ratio, udr2)
}

p8 <- sum(udr2>H3K4m2_dn_ratio)/length(udr2)

d20 <- data.frame(group="Random",ratio=udr2,updn="Near",histone="H3K4me2_M6")

##################################

df <- rbind(d14,d16,d18,d20)
df2 <- rbind(d13,d15,d17,d19)

pdf("NannoEpi_Fig4m2_M4M6_shift.pdf", useDingbats=FALSE, width=3, height=4, family="Arial")

ggboxplot(df, x="histone", y="ratio", fill="updn", palette="jco", size=0.2, bxp.errorbar=T, outlier.shape=NA) + coord_cartesian(ylim=c(0,20)) + 
  geom_point(data=df2, aes(shape=updn), fill="red", colour="red", size=2, position=position_jitterdodge(jitter.width=0)) +
  scale_shape_manual(values=c(24,25)) +
  theme_bw() + theme(text=element_text(size=14, family = "Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("Ratio")

dev.off()

