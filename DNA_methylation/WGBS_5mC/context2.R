#!/usr/bin/env Rscript

# compare genome wide 5mc level among CG, CHG and CHH context

#setwd(dirname(sys.frame(1)$ofile))
setwd('./')

df <- read.delim('5mc_gw.txt',row.names = 1)

df$Total<-NULL

#par(family="Arial",ps=20)

library(ggpubr)
library(reshape2)

library(extrafont)
loadfonts(device = "pdf")

data <- melt(df)

str(data)

pdf("Context5Cgw_fix2.pdf", width=6, height=6, family="Arial", useDingbats = FALSE)

my_comparisons <- list(c("CG","CHG"),c("CHG","CHH"),c("CG","CHH"))

stat.test <- compare_means(value ~ Group, data=data, group.by="variable", method = "wilcox.test")

ggboxplot(data, x="variable", y="value", bxp.errorbar = T, palette="npg", fill="Group", width=0.5) + theme_bw() + 
  xlab("") + ylab("5mC level (%)") + ylim(0.1, 0.185) +
  stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, 
                     method="wilcox.test",label.y = c(0.17,0.175,0.18)) + 
  theme(text=element_text(size=16,  family="Arial"),legend.position="top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  stat_pvalue_manual(
    stat.test, x = "variable", y.position = 0.1,
    label = "p={p}", position = position_dodge(0.8)
  )

dev.off()

