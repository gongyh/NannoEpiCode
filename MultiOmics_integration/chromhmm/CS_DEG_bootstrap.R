#!/usr/bin/env Rscript



data <- read.delim("H0_H24_CS_DEG.txt", header=T, stringsAsFactors=F)
#str(data)

result_df <- data.frame(up=numeric(),up_mean=numeric(),up_sd=numeric(),
                       dn=numeric(),dn_mean=numeric(),dn_sd=numeric())

iter <- 1000

set.seed(123456)

CSa <- c("1","2","3","4","5")
CSb <- c("1","2","3","4","5")

for (CSi in CSa) {
 for (CSj in CSb) {

gene_cs1 <- data[data$CS==as.integer(CSi),]$Gene
gene_cs1_cs1 <- data[data$Transition==paste0(CSi,"_",CSj),]
cs11_up <- length(gene_cs1_cs1[gene_cs1_cs1$DEG=="Up","Gene"])/length(gene_cs1_cs1$Gene)
cs11_dn <- length(gene_cs1_cs1[gene_cs1_cs1$DEG=="Down","Gene"])/length(gene_cs1_cs1$Gene)

ups <- c()
downs <- c()
for (i in 1:iter) {

  sg <- sample(gene_cs1,length(gene_cs1_cs1$Gene))
  sg_all <- data[data$Gene %in% sg,]
  up <- length(sg_all[sg_all$DEG=="Up","Gene"])/length(sg)
  down <- length(sg_all[sg_all$DEG=="Down","Gene"])/length(sg)
  ups <- c(up, ups)
  downs <- c(down, downs)
}

cs_mean_sd <- data.frame(up=cs11_up,up_mean=mean(ups),up_sd=sd(ups),
                         dn=cs11_dn,dn_mean=mean(downs),dn_sd=sd(downs))
row.names(cs_mean_sd) <- paste0("CS_",CSi,"_",CSj)

result_df <- rbind(result_df, cs_mean_sd)

 }
}

result_df
