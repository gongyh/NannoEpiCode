#!/usr/bin/env Rscript

library(reshape2)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

H0_profile <- read.delim("H0.mRp.clN.plotProfile.tab")[,1:902]
bins <- H0_profile[1,3:dim(H0_profile)[2]]
H0 <- H0_profile[2,3:dim(H0_profile)[2]]

H24_profile <- read.delim("H24.mRp.clN.plotProfile.tab")[,1:902]
H24 <- H24_profile[2,3:dim(H24_profile)[2]]

profiles <- rbind(bins=bins, H0=H0, H24=H24)

df <- as.data.frame(t(profiles))

data <- melt(df, id.vars="bins", measure.vars=c("H0","H24"))

pdf("NannoEpi_Fig3A.pdf", useDingbats = FALSE, width = 5, height = 4, family="Arial")

#par(family="Arial", ps=16)
#par(mar=c(5.1, 4.1, 4.1, 3.1), mgp=c(3, 1, 0), las=0)

ggline(data, x="bins", y="value", group="variable", color="variable", palette="npg", 
       plot_type='l', size=0.5) + theme_bw() +
  geom_vline(xintercept=range(200, 700), color='gray', size=0.2, linetype = "dashed") +
  scale_x_discrete(breaks=c(1,100,200,300,400,500,600,700,800,900), 
                  labels=c('-2Kb','-1Kb','TSS','20%','40%','60%','80%','TES','1Kb','2Kb')) +
  theme(text=element_text(size=14,  family="Arial"), 
        plot.margin=margin(t=5, r=15, b=5, l=5, unit = "pt"),
        legend.position=c(0.45,0.4), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Genomic Region (5' -> 3')") + ylab("")

dev.off()

