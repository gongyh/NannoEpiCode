#!/usr/bin/env Rscript

library(reshape2)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

H24_H3k27ac_profile <- read.delim("H24_H3k27ac.plotProfile.tab")[,1:902]
bins <- H24_H3k27ac_profile[1,3:dim(H24_H3k27ac_profile)[2]]
H24_H3k27ac <- (H24_H3k27ac_profile[2,3:dim(H24_H3k27ac_profile)[2]] + H24_H3k27ac_profile[3,3:dim(H24_H3k27ac_profile)[2]])/2

H24_H3k9ac_profile <- read.delim("H24_H3k9ac.plotProfile.tab")[,1:902]
H24_H3k9ac <- (H24_H3k9ac_profile[2,3:dim(H24_H3k9ac_profile)[2]] + H24_H3k9ac_profile[3,3:dim(H24_H3k9ac_profile)[2]])/2

H24_H3k4m2_profile <- read.delim("H24_H3k4m2.plotProfile.tab")[,1:902]
H24_H3k4m2 <- (H24_H3k4m2_profile[2,3:dim(H24_H3k4m2_profile)[2]] + H24_H3k4m2_profile[3,3:dim(H24_H3k4m2_profile)[2]])/2

H24_H3kcr_profile <- read.delim("H24_H3kcr.plotProfile.tab")[,1:902]
H24_H3kcr <- (H24_H3kcr_profile[2,3:dim(H24_H3kcr_profile)[2]] + H24_H3kcr_profile[3,3:dim(H24_H3kcr_profile)[2]])/2

H24_profiles <- rbind(bins=bins, H24_H3k9ac=H24_H3k9ac, H24_H3k27ac=H24_H3k27ac, H24_H3kcr=H24_H3kcr, H24_H3k4m2=H24_H3k4m2)

H24_df <- as.data.frame(t(H24_profiles))

H24_data <- melt(H24_df, id.vars="bins", measure.vars=c("H24_H3k9ac","H24_H3k27ac","H24_H3kcr","H24_H3k4m2"))

pdf("NannoEpi_Fig1B.pdf", useDingbats = FALSE, width = 5, height = 4, family="Arial")

#par(family="Arial", ps=16)
#par(mar=c(5.1, 4.1, 4.1, 3.1), mgp=c(3, 1, 0), las=0)

ggline(H24_data, x="bins", y="value", group="variable", color="variable", palette="npg", 
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=range(200, 700), color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(1,100,200,300,400,500,600,700,800,900), 
                  labels=c('-2Kb','-1Kb','TSS','20%','40%','60%','80%','TES','1Kb','2Kb')) +
  theme(text=element_text(size=14,  family="Arial"), 
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position=c(0.77,0.8), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Genomic Region (5' -> 3')") + ylab("log2(fold change)")

dev.off()

