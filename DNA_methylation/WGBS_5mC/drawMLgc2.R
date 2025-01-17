#!/usr/bin/env Rscript

#setwd(dirname(sys.frame(1)$ofile))
setwd('./')

df <- read.delim('MLgc.txt',row.names = 2)

df$TSSUP.<-df$TSSUP./df$Genome
df$TESDOWN.<-df$TESDOWN./df$Genome
df$Genes.<-df$Genes./df$Genome
df$Exons.<-df$Exons./df$Genome
df$CDSs.<-df$CDSs./df$Genome
df$TEs.<-df$TEs./df$Genome

df$TSSUP<-NULL
df$TESDOWN<-NULL
df$Genes<-NULL
df$Exons<-NULL
df$CDSs<-NULL
df$TEs<-NULL
df$Genome<-NULL

par(family="Arial",ps=20)

library(ggpubr)
library(reshape2)

library(extrafont)
loadfonts(device = "pdf")

data <- melt(df)

stat.test <- compare_means(value ~ Group, data=data, group.by="variable", method = "wilcox.test")

pdf("MLgc_cmp.pdf", width=6, height=6, family="Arial", useDingbats = FALSE)

ggboxplot(data, x="variable", y="value", fill="Group", #position=position_dodge(width=0.8),
          bxp.errorbar = T, palette="npg", width=0.5) + theme_bw() +
  xlab("") + ylab("normalized 5mC levels") + 
  geom_hline(yintercept = 1, linetype="dashed", color = "grey") +
  scale_y_continuous(limits=c(0.95,1.1)) +
  theme(text=element_text(size=16,  family="Arial"),legend.position="top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  stat_pvalue_manual(
    stat.test, x = "variable", y.position = 1.1,
    label = "p={p}", position = position_dodge(0.8)
  )

dev.off()
