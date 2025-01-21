#!/usr/bin/env Rscript

library(reshape2)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

M4_R1_profile <- read.delim("M4_R1.plotProfile.tab")[,1:902]
bins <- M4_R1_profile[1,3:dim(M4_R1_profile)[2]]
M4_R1 <- M4_R1_profile[2,3:dim(M4_R1_profile)[2]]

M4_R2_profile <- read.delim("M4_R2.plotProfile.tab")[,1:902]
M4_R2 <- M4_R2_profile[2,3:dim(M4_R2_profile)[2]]

M4_R3_profile <- read.delim("M4_R3.plotProfile.tab")[,1:902]
M4_R3 <- M4_R3_profile[2,3:dim(M4_R3_profile)[2]]

M6_R1_profile <- read.delim("M6_R1.plotProfile.tab")[,1:902]
M6_R1 <- M6_R1_profile[2,3:dim(M6_R1_profile)[2]]

M6_R2_profile <- read.delim("M6_R2.plotProfile.tab")[,1:902]
M6_R2 <- M6_R2_profile[2,3:dim(M6_R2_profile)[2]]

M6_R3_profile <- read.delim("M6_R3.plotProfile.tab")[,1:902]
M6_R3 <- M6_R3_profile[2,3:dim(M6_R3_profile)[2]]

WT_R1_profile <- read.delim("WT_R1.plotProfile.tab")[,1:902]
WT_R1 <- WT_R1_profile[2,3:dim(WT_R1_profile)[2]]

WT_R2_profile <- read.delim("WT_R2.plotProfile.tab")[,1:902]
WT_R2 <- WT_R2_profile[2,3:dim(WT_R2_profile)[2]]

WT_R3_profile <- read.delim("WT_R3.plotProfile.tab")[,1:902]
WT_R3 <- WT_R3_profile[2,3:dim(WT_R3_profile)[2]]

profiles <- rbind(bins=bins, WT_R1=WT_R1, WT_R2=WT_R2, WT_R3=WT_R3, M4_R1=M4_R1, 
                  M4_R2=M4_R2, M4_R3=M4_R3, M6_R1=M6_R1, M6_R2=M6_R2, M6_R3=M6_R3)

df <- as.data.frame(t(profiles))

data <- melt(df, id.vars="bins", measure.vars=c("WT_R1","WT_R2","WT_R3",
             "M4_R1","M4_R2","M4_R3","M6_R1","M6_R2","M6_R3"))

pdf("NannoEpi_MutProfile.pdf", useDingbats = FALSE, width = 6, height = 4, family="Arial")

colors <- c("gray","gray","gray","red","red","red","blue","blue","blue")

ggline(data, x="bins", y="value", group="variable", color="variable", palette=colors, 
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=range(200, 700), color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(1,100,200,300,400,500,600,700,800,900), 
                  labels=c('-2Kb','-1Kb','TSS','20%','40%','60%','80%','TES','1Kb','2Kb')) +
  theme(text=element_text(size=14,  family="Arial"), 
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position="right", legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Genomic Region (5' -> 3')") + ylab("log2(fold change)")

dev.off()

