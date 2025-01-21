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
H3K9ac_up_ratio <- 162/54

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
H3K9ac_dn_ratio <- 239/168

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
H3Kcr_up_ratio <- 312/79

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
H3Kcr_dn_ratio <- 454/141

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
H3K27ac_up_ratio <- 299/56

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
H3K27ac_dn_ratio <- 461/123

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
H3K4m2_peaks <- read.delim("H3k4m2/H3k4m2_DEPeaks.txt", header=T, stringsAsFactors=F)
num <- 145
H3K4m2_up_ratio <- 12/68

d13 <- data.frame(group="Peaks_DE",ratio=H3K4m2_up_ratio,updn="Up",histone="H3K4me2")

udr <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$Gene, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(ratio, udr)
}

p7 <- sum(udr>H3K4m2_up_ratio)/length(udr)

d14 <- data.frame(group="Random",ratio=udr,updn="Up",histone="H3K4me2")

num <- 118
H3K4m2_dn_ratio <- 9/58

d15 <- data.frame(group="Peaks_DE",ratio=H3K4m2_dn_ratio,updn="Down",histone="H3K4me2")

udr2 <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$Gene, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg_all[sg_all$DEG=="Up","Gene"])
  udr2 <- c(ratio, udr2)
}

p8 <- sum(udr2>H3K4m2_dn_ratio)/length(udr2)

d16 <- data.frame(group="Random",ratio=udr2,updn="Down",histone="H3K4me2")

##################################

df <- rbind(d2,d4,d6,d8,d10,d12,d14,d16)
df2 <- rbind(d1,d3,d5,d7,d9,d11,d13,d15)

pdf("NannoEpi_Fig3C.pdf", useDingbats=FALSE, width=8, height=4, family="Arial")

ggbarplot(df, x="histone", y="ratio", fill="updn", palette="jco", width=0.5, alpha=0.5,
          add="mean_sd", group="updn", color="updn", position=position_dodge(0.7)) + 
  geom_point(data=df2, aes(shape=updn), color="red", size=2,
             position=position_jitterdodge(jitter.width=0, dodge.width=0.7)) +
  theme_bw() + theme(text=element_text(size=14, family = "Arial"), legend.position=c(0.92,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("Ratio")

dev.off()

