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
H3K4m2_peaks <- read.delim("../H24vsH0_H3K4me2_shift.txt", header=F, stringsAsFactors=F)
num <- 1274
H3K4m2_up_ratio <- 398/121

d13 <- data.frame(group="Peaks_DE",ratio=H3K4m2_up_ratio,updn="Far",histone="H3K4me2")

udr <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$V1, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr <- c(ratio, udr)
}

p7 <- sum(udr>H3K4m2_up_ratio)/length(udr)

d14 <- data.frame(group="Random",ratio=udr,updn="Far",histone="H3K4me2")

num <- 790
H3K4m2_dn_ratio <- 93/243

d15 <- data.frame(group="Peaks_DE",ratio=H3K4m2_dn_ratio,updn="Near",histone="H3K4me2")

udr2 <- c()
for (i in 1:iter) {

  sg <- sample(H3K4m2_peaks$V1, num)
  sg_all <- data[data$Gene %in% sg,]
  ratio <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg_all[sg_all$DEG=="Down","Gene"])
  udr2 <- c(ratio, udr2)
}

p8 <- sum(udr2>H3K4m2_dn_ratio)/length(udr2)

d16 <- data.frame(group="Random",ratio=udr2,updn="Near",histone="H3K4me2")

##################################

df <- rbind(d14,d16)
df2 <- rbind(d13,d15)

df$ratio2 <- log2(df$ratio)
df2$ratio2 <- log2(df2$ratio)

pdf("NannoEpi_Fig4me2.pdf", useDingbats=FALSE, width=3, height=3, family="Arial")

ggboxplot(df, x="histone", y="ratio2", fill="updn", palette="jco", size=0.2, bxp.errorbar=T, facet.by="updn", nrow=1) + 
  geom_point(data=df2, aes(shape=updn), fill="red", colour="red", size=2, position=position_jitterdodge(jitter.width=0)) +
  scale_shape_manual(values=c(24,25)) +
  theme_bw() + theme(text=element_text(size=14, family = "Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("log2 ratio (Up/Down)")

dev.off()

