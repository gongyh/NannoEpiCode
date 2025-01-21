#!/usr/bin/env Rscript

library(ggpubr)
library(reshape2)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

H0 <- read.table('H0_phases.txt')
H24 <- read.table('H24_phases.txt')

H0_hist<-hist(H0[H0$V1<=300,],breaks=1:300,plot=FALSE)
H24_hist<-hist(H24[H24$V1<=300,],breaks=1:300,plot=FALSE)

phases_df <- data.frame(Phase=100:299,H0=H0_hist$density[100:299],H24=H24_hist$density[100:299])

phases <- melt(phases_df, id.vars="Phase")

pdf("NannoEpi_Fig3_phases.pdf", useDingbats = FALSE, width = 4, height = 3, family="Arial")

ggline(phases,x="Phase",y="value",color="variable",palette="npg", plot_type="l") + theme_bw() +
  geom_vline(xintercept=178, color='gray',linetype='dashed', size=0.2) +
  scale_x_continuous(breaks=c(100,150,178,200,250,300),labels=c(100,150,178,200,250,300)) + 
  theme(text=element_text(size=11, family="Arial"), 
        legend.position=c(0.8,0.8),legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Phase (bp)") + ylab("Phase density")

dev.off()

