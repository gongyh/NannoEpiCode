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
H3K4m2_peaks <- read.delim("H3k4m2/M4_DEPeaks.txt", header=T, stringsAsFactors=F)
num <- 1536
H3K4m2_up_ratio <- 163/16

d13 <- data.frame(group="Peaks_DE",ratio=1/H3K4m2_up_ratio,updn="Up",histone="H3K4me2_M4")

udr <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$Gene, num)
  sg_all <- data4[data4$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(1/ratio, udr)
}

p7 <- sum(udr>H3K4m2_up_ratio)/length(udr)
p7

d14 <- data.frame(group="Random",ratio=udr,updn="Up",histone="H3K4me2_M4")

num <- 1274
H3K4m2_dn_ratio <- 10/42

d15 <- data.frame(group="Peaks_DE",ratio=1/H3K4m2_dn_ratio,updn="Down",histone="H3K4me2_M4")

udr2 <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$Gene, num)
  sg_all <- data4[data4$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg_all[sg_all$DEG=="Up","Gene"])
  udr2 <- c(1/ratio, udr2)
}

p8 <- sum(udr2>H3K4m2_dn_ratio)/length(udr2)
p8

d16 <- data.frame(group="Random",ratio=udr2,updn="Down",histone="H3K4me2_M4")

##################################
H3K4m2_peaks <- read.delim("H3k4m2/M6_DEPeaks.txt", header=T, stringsAsFactors=F)
num <- 1906
H3K4m2_up_ratio <- 93/26

d17 <- data.frame(group="Peaks_DE",ratio=1/H3K4m2_up_ratio,updn="Up",histone="H3K4me2_M6")

udr <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$Gene, num)
  sg_all <- data6[data6$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(1/ratio, udr)
}

p7 <- sum(udr>H3K4m2_up_ratio)/length(udr)
p7

d18 <- data.frame(group="Random",ratio=udr,updn="Up",histone="H3K4me2_M6")

num <- 1733
H3K4m2_dn_ratio <- 37/50

d19 <- data.frame(group="Peaks_DE",ratio=1/H3K4m2_dn_ratio,updn="Down",histone="H3K4me2_M6")

udr2 <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$Gene, num)
  sg_all <- data6[data6$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg_all[sg_all$DEG=="Up","Gene"])
  udr2 <- c(1/ratio, udr2)
}

p8 <- sum(udr2>H3K4m2_dn_ratio)/length(udr2)
p8

d20 <- data.frame(group="Random",ratio=udr2,updn="Down",histone="H3K4me2_M6")

##################################

df <- rbind(d14, d16, d18, d20)
df2 <- rbind(d13, d15, d17, d19)

pdf("NannoEpi_M4M6_2.pdf", useDingbats=FALSE, width=4, height=4, family="Arial")

ggboxplot(df, x="histone", y="ratio", fill="updn", palette="jco", size=0.2, bxp.errorbar=T, outlier.shape=NA) + coord_cartesian(ylim=c(0,20)) +
  geom_point(data=df2, aes(shape=updn), fill="red", colour="red", size=2, position=position_jitterdodge(jitter.width=0)) +
  scale_shape_manual(values=c(24,25)) +
  theme_bw() + theme(text=element_text(size=14, family = "Arial"), legend.position=c(0.85,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("1/Ratio")

dev.off()

