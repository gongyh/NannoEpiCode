#!/usr/bin/env Rscript

library(stats)

data <- read.delim("H0_H24_CS_DEG.txt", header=T, stringsAsFactors=F)
#str(data)

result_df <- data.frame(raw=character(),cs=character(),deg=character(),
                        realp=numeric(),overalp=numeric(),pvalue=numeric())

set.seed(123456)

CSa <- c("1","2","3","4","5")
CSb <- c("1","2","3","4","5")

for (CSi in CSa) {
 for (CSj in CSb) {

cs1 <- data[data$CS==as.integer(CSi),]
p_up <- length(cs1[cs1$DEG=="Up","Gene"])/length(cs1$Gene)
p_dn <- length(cs1[cs1$DEG=="Down","Gene"])/length(cs1$Gene)
gene_cs1 <- cs1$gene
gene_cs1_cs1 <- data[data$Transition==paste0(CSi,"_",CSj),]
cs11 <- length(gene_cs1_cs1$Gene)
cs11_up <- length(gene_cs1_cs1[gene_cs1_cs1$DEG=="Up","Gene"])
cs11_dn <- length(gene_cs1_cs1[gene_cs1_cs1$DEG=="Down","Gene"])


if (cs11>0) {
bt_up <- binom.test(x=cs11_up,n=cs11,p=p_up)
bt_dn <- binom.test(x=cs11_dn,n=cs11,p=p_dn)
} else {
bt_up$p.value=NA
bt_dn$p.value=NA
}

up_test <- data.frame(raw=CSi,cs=paste0("->",CSj),deg="Up",realp=100*cs11_up/cs11,overalp=100*p_up,pvalue=bt_up$p.value)
dn_test <- data.frame(raw=CSi,cs=paste0("->",CSj),deg="Down",realp=100*cs11_dn/cs11,overalp=100*p_dn,pvalue=bt_dn$p.value)

#row.names(cs_mean_sd) <- paste0("CS_",CSi,"_",CSj)

result_df <- rbind(result_df, up_test, dn_test)

 }
}


result_df
#q()
total_df <- result_df[result_df$cs=="->1",c("raw","deg","overalp")]

library(ggpubr)
library(reshape2)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

pdf("NannoEpi_Fig5_binom.pdf", useDingbats = FALSE, width = 4, height = 5, family="Arial")

ggbarplot(result_df,x="cs",y="realp",facet.by=c("raw","deg"),scales="free",color="white",fill="cs",palette="npg") + 
  theme_bw() + scale_y_continuous(expand = expansion(mult=c(0,.1))) +
#  geom_hline(data=total_df, aes(yintercept="overalp"), color='gray',linetype='dashed', size=0.2) +
  theme(text=element_text(size=11, family="Arial"), 
        legend.position="none",legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Chromatin state transition") + ylab("Percent of DEG")

dev.off()

