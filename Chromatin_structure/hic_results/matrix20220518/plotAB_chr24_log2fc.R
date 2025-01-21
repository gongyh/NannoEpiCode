#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

HC_R1_df <- read.table('NannoH0_compartments.cis.E1.bedgraph', sep="\t", header=F)
#HC_R2_df <- read.table('H0rep_compartments.cis.E1.bedgraph', sep="\t", header=F)
LC_R1_df <- read.table('NannoH24_compartments.cis.E1.bedgraph', sep="\t", header=F)
#LC_R2_df <- read.table('H24rep_compartments.cis.E1.bedgraph', sep="\t", header=F)

bin_gene <- read.table('bins_10000_genes.bed', sep="\t", header=F)
gene_log2fc <- read.table('LC_logfc.tsv', sep="\t", header=F)

gene_tpm <- read.table('HC_VLC.TMM.TPM.matrix.noCpMt',sep="\t", header=T)

bin_gene_tpm <- merge(bin_gene,gene_log2fc,by.x="V4",by.y="V1")
#str(bin_gene_log2fc)
bin_gene_log2fc <- merge(bin_gene_tpm, gene_tpm, by.x="V4", by.y="Genes")
bin_gene_log2fc$pos <- (bin_gene_log2fc$V2.x + bin_gene_log2fc$V3)/2

HC_R1_df$pos <- (HC_R1_df$V2 + HC_R1_df$V3)/2
#HC_R2_df$pos <- (HC_R2_df$V2 + HC_R2_df$V3)/2
LC_R1_df$pos <- (LC_R1_df$V2 + LC_R1_df$V3)/2
#LC_R2_df$pos <- (LC_R2_df$V2 + LC_R2_df$V3)/2

#HC_R1_df$V4 <- 0-HC_R1_df$V4
#LC_R1_df$V4 <- 0-LC_R1_df$V4

HC_R1_df$sign <- HC_R1_df$V4 >= 0.0
#HC_R2_df$sign <- HC_R2_df$V4 >= 0.0
LC_R1_df$sign <- LC_R1_df$V4 >= 0.0
#LC_R2_df$sign <- LC_R2_df$V4 >= 0.0

pdf("NannoEpi_HiC_AB_chr24_log2fc.pdf", useDingbats=FALSE, width=5, height=5, family="Arial", onefile=F)

HC_data <- HC_R1_df[HC_R1_df$V1=="chr24",]
LC_data <- LC_R1_df[LC_R1_df$V1=="chr24",]

bgl_data <- bin_gene_log2fc[bin_gene_log2fc$V1=="chr24",]
#str(bgl_data)

HC_data$grp <- 'HC'
LC_data$grp <- 'LC'

data <- rbind(HC_data, LC_data)

a <- ggdotchart(data, x="pos", y="V4", sorting="none", add = "segment", add.params=list(size=0.5,color="sign"),
           facet.by="grp", ncol=1, strip.position='right', size=0.5, color="sign", dot.size=0.25, palette=c("blue","red")) +
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black")) +
  xlab("Chromosome 24") + ylab("Eigenvector")


bgl_data$HC <- bgl_data$Control
bgl_data$LC <- bgl_data$VLC
ab24 <- melt(bgl_data,measure.vars=c("HC","LC"),id.vars=c("V1","pos"),variable.name="grp")
ab24_df <- merge(ab24,data)
#ab24_df$value <- ab24_df$value+1
ab24_df$x <- paste(ab24_df$grp,ab24_df$sign,sep="_")
ab24_df$x <- factor(ab24_df$x,levels=c("HC_TRUE","HC_FALSE","LC_TRUE","LC_FALSE"),ordered=T)
#str(ab24_df)

cmps<-list(c("HC_TRUE","HC_FALSE"),c("LC_TRUE","LC_FALSE"))

ab <- ggbarplot(ab24_df, x="x", y="value", group="sign", fill="sign", palette=c("blue","red"), 
                add="mean_se", error.plot="errorbar", ylim=c(0,425)) +
#  geom_hline(yintercept=0, linetype="dashed", colour="gray", size=0.01) + 
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("RNA-Seq (TPM)") + 
  stat_compare_means(comparisons=cmps, method.args=list(alternative="greater"), method="t.test", 
                     label.y=c(300,410), tip.length=0.003)

data2 <- merge(HC_data,LC_data,by="pos")
df <- data2[,c("pos","sign.x","sign.y")]
df$type <- paste(data2$sign.x,data2$sign.y,sep="_")
data3 <- merge(bgl_data,df)
data3$type <- factor(data3$type, levels=c("TRUE_TRUE","TRUE_FALSE","FALSE_TRUE","FALSE_FALSE"), ordered=T)
cmps2 <- list(c("TRUE_FALSE","FALSE_TRUE"))

b <- ggboxplot(data3, x="type", y="V2.y", fill="type", pallette="jco") +
#  geom_hline(yintercept=0, linetype="dashed", colour="gray", size=0.01) + 
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("RNA-Seq (log2fc)") +
  stat_compare_means(comparisons=cmps2, method="t.test", label.y=6)

ggarrange(a,ggarrange(ab,b,nrow=1),ncol=1)

dev.off()

