#!/usr/bin/env Rscript

library(reshape2)
library(ggpubr)
library(stringr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

H0_H3k27ac_profile <- read.delim("../H0_H3k27ac.plotProfile.tab")[,1:902]
bins1 <- H0_H3k27ac_profile[1,3:dim(H0_H3k27ac_profile)[2]]
H0_H3k27ac_R1 <- H0_H3k27ac_profile[2,3:dim(H0_H3k27ac_profile)[2]]
H0_H3k27ac_R2 <- H0_H3k27ac_profile[3,3:dim(H0_H3k27ac_profile)[2]]

H0_H3k9ac_profile <- read.delim("../H0_H3k9ac.plotProfile.tab")[,1:902]
H0_H3k9ac_R1 <- H0_H3k9ac_profile[2,3:dim(H0_H3k9ac_profile)[2]]
H0_H3k9ac_R2 <- H0_H3k9ac_profile[3,3:dim(H0_H3k9ac_profile)[2]]

H0_H3k4m2_profile <- read.delim("../H0_H3k4m2.plotProfile.tab")[,1:902]
H0_H3k4m2_R1 <- H0_H3k4m2_profile[2,3:dim(H0_H3k4m2_profile)[2]]
H0_H3k4m2_R2 <- H0_H3k4m2_profile[3,3:dim(H0_H3k4m2_profile)[2]]

H0_H3kcr_profile <- read.delim("../H0_H3kcr.plotProfile.tab")[,1:902]
H0_H3kcr_R1 <- H0_H3kcr_profile[2,3:dim(H0_H3kcr_profile)[2]]
H0_H3kcr_R2 <- H0_H3kcr_profile[3,3:dim(H0_H3kcr_profile)[2]]

H0_profiles <- rbind(bins=bins1, H0_H3k9ac=H0_H3k9ac_R1, H0_H3k27ac=H0_H3k27ac_R1, H0_Kcr=H0_H3kcr_R1, H0_H3k4me2=H0_H3k4m2_R1)

H24_H3k27ac_profile <- read.delim("../H24_H3k27ac.plotProfile.tab")[,1:902]
bins2 <- H24_H3k27ac_profile[1,3:dim(H24_H3k27ac_profile)[2]]
H24_H3k27ac_R1 <- H24_H3k27ac_profile[2,3:dim(H24_H3k27ac_profile)[2]]
H24_H3k27ac_R2 <- H24_H3k27ac_profile[3,3:dim(H24_H3k27ac_profile)[2]]

H24_H3k9ac_profile <- read.delim("../H24_H3k9ac.plotProfile.tab")[,1:902]
H24_H3k9ac_R1 <- H24_H3k9ac_profile[2,3:dim(H24_H3k9ac_profile)[2]] 
H24_H3k9ac_R2 <- H24_H3k9ac_profile[3,3:dim(H24_H3k9ac_profile)[2]]

H24_H3k4m2_profile <- read.delim("../H24_H3k4m2.plotProfile.tab")[,1:902]
H24_H3k4m2_R1 <- H24_H3k4m2_profile[2,3:dim(H24_H3k4m2_profile)[2]]
H24_H3k4m2_R2 <- H24_H3k4m2_profile[3,3:dim(H24_H3k4m2_profile)[2]]

H24_H3kcr_profile <- read.delim("../H24_H3kcr.plotProfile.tab")[,1:902]
H24_H3kcr_R1 <- H24_H3kcr_profile[2,3:dim(H24_H3kcr_profile)[2]]
H24_H3kcr_R2 <- H24_H3kcr_profile[3,3:dim(H24_H3kcr_profile)[2]]

H24_profiles <- rbind(bins=bins2, H24_H3k9ac=H24_H3k9ac_R1, H24_H3k27ac=H24_H3k27ac_R1, H24_Kcr=H24_H3kcr_R1, H24_H3k4me2=H24_H3k4m2_R1)

H0_df <- as.data.frame(t(H0_profiles))
H24_df <- as.data.frame(t(H24_profiles))

df <- merge(H0_df,H24_df,by="bins")

#str(df)

R1_data <- melt(df, id.vars="bins", measure.vars=c("H0_H3k9ac","H0_H3k27ac","H0_Kcr","H0_H3k4me2","H24_H3k9ac","H24_H3k27ac","H24_Kcr","H24_H3k4me2"))

out <- str_split_fixed(R1_data$variable, "_", 2)
R1_data$condition <- out[,1]
R1_data$histone <- out[,2]

str(R1_data)

pdf("NannoEpi_histone_profile_R1.pdf", useDingbats = FALSE, width = 6, height = 6, family="Arial")

#par(family="Arial", ps=16)
#par(mar=c(5.1, 4.1, 4.1, 3.1), mgp=c(3, 1, 0), las=0)

R1_data$histone <- factor(R1_data$histone, levels=c("H3k9ac","H3k27ac","Kcr","H3k4me2"), ordered=T)

ggline(R1_data, x="bins", y="value", group="condition", color="condition", palette=c("#666666","red"), facet.by="histone", nrow=2,
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=range(200, 700), color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(50,200,450,700,850), 
                  labels=c('-1.5K','TSS','50%','TES','1.5K')) +
  theme(text=element_text(size=14,  family="Arial"), 
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position=c(0.25,0.9), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Genomic Region (5' -> 3')") + ylab("log2(fold change)")

dev.off()

