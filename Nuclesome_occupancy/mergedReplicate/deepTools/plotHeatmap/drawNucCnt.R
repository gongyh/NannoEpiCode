#!/usr/bin/env Rscript

library(reshape2)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

HC_nuc_cnt <- read.delim("HC_nucleosomes.txt")
HC_nuc_cnt$pos <- rownames(HC_nuc_cnt)
data <- melt(HC_nuc_cnt)
str(data)

LC_nuc_cnt <- read.delim("LC_nucleosomes.txt")
LC_nuc_cnt$pos <- rownames(LC_nuc_cnt)
data2 <- melt(LC_nuc_cnt)
str(data2)

HC_nuc_cnt2 <- read.delim("HC_nucleosomes_TES.txt")
HC_nuc_cnt2$pos <- rownames(HC_nuc_cnt2)
data11 <- melt(HC_nuc_cnt2)
str(data11)

LC_nuc_cnt2 <- read.delim("LC_nucleosomes_TES.txt")
LC_nuc_cnt2$pos <- rownames(LC_nuc_cnt2)
data22 <- melt(LC_nuc_cnt2)
str(data22)

HC_nuc_cnt2k <- read.delim("HC_nucleosomes2.txt")
HC_nuc_cnt2k$pos <- rownames(HC_nuc_cnt2k)
data2k <- melt(HC_nuc_cnt2k)
str(data2k)

LC_nuc_cnt2k <- read.delim("LC_nucleosomes2.txt")
LC_nuc_cnt2k$pos <- rownames(LC_nuc_cnt2k)
data2kL <- melt(LC_nuc_cnt2k)
str(data2kL)

pdf("NannoEpi_NucCnt.pdf", useDingbats = FALSE, width = 5, height = 4, family="Arial")

ggline(data, x="pos", y="value", group="variable", color="variable", palette="npg", 
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=200, color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(1,100,200,300,400), labels=c('-2Kb','-1Kb','TSS','1Kb','2Kb')) +
  theme(text=element_text(size=14,  family="Arial"), 
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position=c(0.4,0.82), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("# nucleosomes (1kb bins)")

ggline(data2, x="pos", y="value", group="variable", color="variable", palette="npg",
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=200, color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(1,100,200,300,400), labels=c('-2Kb','1Kb','TSS','1Kb','2Kb')) +
  theme(text=element_text(size=14,  family="Arial"),
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position=c(0.4,0.82), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("# nucleosomes (1kb bins)")

ggline(data11, x="pos", y="value", group="variable", color="variable", palette="npg",
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=200, color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(1,100,200,300,400), labels=c('-2Kb','1Kb','TES','1Kb','2Kb')) +
  theme(text=element_text(size=14,  family="Arial"),
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position=c(0.55,0.82), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("# nucleosomes (1kb bins)")

ggline(data22, x="pos", y="value", group="variable", color="variable", palette="npg",
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=200, color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(1,100,200,300,400), labels=c('-2Kb','1Kb','TES','1Kb','2Kb')) +
  theme(text=element_text(size=14,  family="Arial"),
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position=c(0.55,0.82), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("# nucleosomes (1kb bins)")

ggline(data2k, x="pos", y="value", group="variable", color="variable", palette="npg",
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=200, color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(1,100,200,300,400), labels=c('-2Kb','-1Kb','TSS','1Kb','2Kb')) +
  theme(text=element_text(size=14,  family="Arial"),
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position=c(0.4,0.82), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("# nucleosomes (2kb bins)")

ggline(data2kL, x="pos", y="value", group="variable", color="variable", palette="npg",
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=200, color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(1,100,200,300,400), labels=c('-2Kb','-1Kb','TSS','1Kb','2Kb')) +
  theme(text=element_text(size=14,  family="Arial"),
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position=c(0.45,0.18), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("# nucleosomes (2kb bins)")

dev.off()

