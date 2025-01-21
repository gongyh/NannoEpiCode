#!/usr/bin/env Rscript

setwd('./')

set.seed(123456)
iter <- 1000

library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

data <- read.delim("H0_H24_DEG.txt", header=T, stringsAsFactors=F)

##################################
H3K9ac_peaks <- read.delim("H3k9ac/H3k9ac_DEPeaks.txt", header=T, stringsAsFactors=F)
num <- 774
H3K9ac_up_ratio <- 162/168

d1 <- data.frame(group="Peaks_DE",ratio=H3K9ac_up_ratio,updn="Up",histone="H3K9ac")

udr <- c()
for (i in 1:iter) {
  
  sg <- sample(H3K9ac_peaks$Gene, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(ratio, udr)
}

p1 <- sum(udr>H3K9ac_up_ratio)/length(udr)

d2 <- data.frame(group="Random",ratio=udr,updn="Up",histone="H3K9ac")

num <- 1049
H3K9ac_dn_ratio <- 239/54

d3 <- data.frame(group="Peaks_DE",ratio=H3K9ac_dn_ratio,updn="Down",histone="H3K9ac")

udr2 <- c()
for (i in 1:iter) {
  
  sg <- sample(H3K9ac_peaks$Gene, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg_all[sg_all$DEG=="Up","Gene"])
  udr2 <- c(ratio, udr2)
}

p2 <- sum(udr2>H3K9ac_dn_ratio)/length(udr2)

d4 <- data.frame(group="Random",ratio=udr2,updn="Down",histone="H3K9ac")

##################################
H3Kcr_peaks <- read.delim("H3kcr/H3kcr_DEPeaks.txt", header=T, stringsAsFactors=F)
num <- 1163
H3Kcr_up_ratio <- 312/141

d9 <- data.frame(group="Peaks_DE",ratio=H3Kcr_up_ratio,updn="Up",histone="H3Kcr")

udr <- c()
for (i in 1:iter) {
  
  sg <- sample(H3Kcr_peaks$Gene, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(ratio, udr)
}

p5 <- sum(udr>H3Kcr_up_ratio)/length(udr)

d10 <- data.frame(group="Random",ratio=udr,updn="Up",histone="H3Kcr")

num <- 1289
H3Kcr_dn_ratio <- 454/79

d11 <- data.frame(group="Peaks_DE",ratio=H3Kcr_dn_ratio,updn="Down",histone="H3Kcr")

udr2 <- c()
for (i in 1:iter) {
  
  sg <- sample(H3Kcr_peaks$Gene, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg_all[sg_all$DEG=="Up","Gene"])
  udr2 <- c(ratio, udr2)
}

p6 <- sum(udr2>H3Kcr_dn_ratio)/length(udr2)

d12 <- data.frame(group="Random",ratio=udr2,updn="Down",histone="H3Kcr")

##################################
H3K27ac_peaks <- read.delim("H3k27ac/H3k27ac_DEPeaks.txt", header=T, stringsAsFactors=F)
num <- 958
H3K27ac_up_ratio <- 299/123

d5 <- data.frame(group="Peaks_DE",ratio=H3K27ac_up_ratio,updn="Up",histone="H3K27ac")

udr <- c()
for (i in 1:iter) {
  
  sg <- sample(H3K27ac_peaks$Gene, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(ratio, udr)
}

p3 <- sum(udr>H3K27ac_up_ratio)/length(udr)

d6 <- data.frame(group="Random",ratio=udr,updn="Up",histone="H3K27ac")

num <- 1211
H3K27ac_dn_ratio <- 461/56

d7 <- data.frame(group="Peaks_DE",ratio=H3K27ac_dn_ratio,updn="Down",histone="H3K27ac")

udr2 <- c()
for (i in 1:iter) {
  
  sg <- sample(H3K27ac_peaks$Gene, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg_all[sg_all$DEG=="Up","Gene"])
  udr2 <- c(ratio, udr2)
}

p4 <- sum(udr2>H3K27ac_dn_ratio)/length(udr2)

d8 <- data.frame(group="Random",ratio=udr2,updn="Down",histone="H3K27ac")

##################################
H3K4m2_peaks <- read.delim("../H24vsH0_H3K4me2_shift.txt", header=F, stringsAsFactors=F)
num <- 1274
H3K4m2_up_ratio <- 398/121

d13 <- data.frame(group="Peaks_DE",ratio=H3K4m2_up_ratio,updn="Up",histone="H3K4me2")

udr <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$V1, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(ratio, udr)
}

p7 <- sum(udr>H3K4m2_up_ratio)/length(udr)
p7
d14 <- data.frame(group="Random",ratio=udr,updn="Up",histone="H3K4me2")

num <- 790
H3K4m2_dn_ratio <- 243/93

d15 <- data.frame(group="Peaks_DE",ratio=H3K4m2_dn_ratio,updn="Down",histone="H3K4me2")

udr2 <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$V1, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg_all[sg_all$DEG=="Up","Gene"])
  udr2 <- c(ratio, udr2)
}

p8 <- sum(udr2>H3K4m2_dn_ratio)/length(udr2)
p8
d16 <- data.frame(group="Random",ratio=udr2,updn="Down",histone="H3K4me2")

##################################

df <- rbind(d1,d1,d1,d2,d3,d3,d3,d4,d5,d5,d5,d6,d7,d7,d7,d8,d9,d9,d9,d10,
            d11,d11,d11,d12,d13,d13,d13,d14,d15,d15,d15,d16)

pdf("NannoEpi_Fig6.pdf", useDingbats=FALSE, width=8, height=4, family="Arial")

ggbarplot(df, x="group", y="ratio", fill="lightgray", palette="npg", width=0.5, facet.by="histone",
          add="mean_sd", group="updn", color="updn", position=position_dodge(0.7), nrow=1) + 
  guides(color=guide_legend(title="Peaks")) +
  theme_bw() + theme(text=element_text(size=14, family = "Arial"), legend.position="right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("Ratio")

dev.off()

