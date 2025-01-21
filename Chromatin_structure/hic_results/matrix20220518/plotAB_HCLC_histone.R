#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

HC_df <- read.table('HC_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)
LC_df <- read.table('LC_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)

colnames(HC_df) <- c("chr","start","end","E1")
colnames(LC_df) <- c("chr","start","end","E1")

HC_df$pos <- (HC_df$start + HC_df$end)/2
LC_df$pos <- (LC_df$start + LC_df$end)/2

HC_df$sign <- HC_df$E1 >= 0.0
LC_df$sign <- LC_df$E1 >= 0.0

HC_df$grp <- 'HC'
LC_df$grp <- 'LC'

data1 <- merge(HC_df,LC_df,by=c("chr","start","end","pos"))

cmps2 <- list(c("TRUE_FALSE","TRUE_TRUE"),c("FALSE_TRUE","FALSE_FALSE"))

chosed <- c("chr1","chr2","chr3","chr4","chr5","chr9","chr10",
            "chr11","chr12","chr15","chr16","chr17","chr18","chr19",
            "chr23","chr24","chr25","chr27","chr28","chr30")

# H3k9ac
bin_H3k9ac <- read.table('bins_40000_H3k9ac.bed', sep="\t", header=F)
colnames(bin_H3k9ac) <- c("chr","start","end","interval")
H3k9ac_log2fc <- read.table('H3k9ac.consensus_peaks.results.txt', sep="\t", header=T)[,c(1,8)]
colnames(H3k9ac_log2fc) <- c("interval","log2fc")
bin_H3k9ac <- merge(bin_H3k9ac,H3k9ac_log2fc,by="interval")
bin_H3k9ac$pos <- (bin_H3k9ac$start + bin_H3k9ac$end)/2
data_H3k9ac <- merge(bin_H3k9ac, data1, by=c("chr","start","end","pos"))
data_H3k9ac$chrom <- as.numeric(sub("chr","",data_H3k9ac$chr))
data_H3k9ac$type <- paste(data_H3k9ac$sign.x,data_H3k9ac$sign.y,sep="_")
data_H3k9ac <- data_H3k9ac[data_H3k9ac$E1.x*data_H3k9ac$E1.y != 0,]
data_H3k9ac$type <- factor(data_H3k9ac$type, levels=c("TRUE_FALSE","TRUE_TRUE","FALSE_TRUE","FALSE_FALSE"), ordered=T)
data_H3k9ac <- data_H3k9ac[data_H3k9ac$chr %in% chosed,]
data_H3k9ac$histone <- 'H3K9ac'

# H3k27ac
bin_H3k27ac <- read.table('bins_40000_H3k27ac.bed', sep="\t", header=F)
colnames(bin_H3k27ac) <- c("chr","start","end","interval")
H3k27ac_log2fc <- read.table('H3k27ac.consensus_peaks.results.txt', sep="\t", header=T)[,c(1,8)]
colnames(H3k27ac_log2fc) <- c("interval","log2fc")
bin_H3k27ac <- merge(bin_H3k27ac,H3k27ac_log2fc,by="interval")
bin_H3k27ac$pos <- (bin_H3k27ac$start + bin_H3k27ac$end)/2
data_H3k27ac <- merge(bin_H3k27ac, data1, by=c("chr","start","end","pos"))
data_H3k27ac$chrom <- as.numeric(sub("chr","",data_H3k27ac$chr))
data_H3k27ac$type <- paste(data_H3k27ac$sign.x,data_H3k27ac$sign.y,sep="_")
data_H3k27ac <- data_H3k27ac[data_H3k27ac$E1.x*data_H3k27ac$E1.y != 0,]
data_H3k27ac$type <- factor(data_H3k27ac$type, levels=c("TRUE_FALSE","TRUE_TRUE","FALSE_TRUE","FALSE_FALSE"), ordered=T)
data_H3k27ac <- data_H3k27ac[data_H3k27ac$chr %in% chosed,]
data_H3k27ac$histone <- 'H3K27ac'

# H3kcr
bin_H3kcr <- read.table('bins_40000_H3kcr.bed', sep="\t", header=F)
colnames(bin_H3kcr) <- c("chr","start","end","interval")
H3kcr_log2fc <- read.table('H3kcr.consensus_peaks.results.txt', sep="\t", header=T)[,c(1,8)]
colnames(H3kcr_log2fc) <- c("interval","log2fc")
bin_H3kcr <- merge(bin_H3kcr,H3kcr_log2fc,by="interval")
bin_H3kcr$pos <- (bin_H3kcr$start + bin_H3kcr$end)/2
data_H3kcr <- merge(bin_H3kcr, data1, by=c("chr","start","end","pos"))
data_H3kcr$chrom <- as.numeric(sub("chr","",data_H3kcr$chr))
data_H3kcr$type <- paste(data_H3kcr$sign.x,data_H3kcr$sign.y,sep="_")
data_H3kcr <- data_H3kcr[data_H3kcr$E1.x*data_H3kcr$E1.y != 0,]
data_H3kcr$type <- factor(data_H3kcr$type, levels=c("TRUE_FALSE","TRUE_TRUE","FALSE_TRUE","FALSE_FALSE"), ordered=T)
data_H3kcr <- data_H3kcr[data_H3kcr$chr %in% chosed,]
data_H3kcr$histone <- 'Kcr'

# H3k4m2
bin_H3k4m2 <- read.table('bins_40000_H3k4m2.bed', sep="\t", header=F)
colnames(bin_H3k4m2) <- c("chr","start","end","interval")
H3k4m2_log2fc <- read.table('H3k4m2.consensus_peaks.results.txt', sep="\t", header=T)[,c(1,8)]
colnames(H3k4m2_log2fc) <- c("interval","log2fc")
bin_H3k4m2 <- merge(bin_H3k4m2,H3k4m2_log2fc,by="interval")
bin_H3k4m2$pos <- (bin_H3k4m2$start + bin_H3k4m2$end)/2
data_H3k4m2 <- merge(bin_H3k4m2, data1, by=c("chr","start","end","pos"))
data_H3k4m2$chrom <- as.numeric(sub("chr","",data_H3k4m2$chr))
data_H3k4m2$type <- paste(data_H3k4m2$sign.x,data_H3k4m2$sign.y,sep="_")
data_H3k4m2 <- data_H3k4m2[data_H3k4m2$E1.x*data_H3k4m2$E1.y != 0,]
data_H3k4m2$type <- factor(data_H3k4m2$type, levels=c("TRUE_FALSE","TRUE_TRUE","FALSE_TRUE","FALSE_FALSE"), ordered=T)
data_H3k4m2 <- data_H3k4m2[data_H3k4m2$chr %in% chosed,]
data_H3k4m2$histone <- 'H3K4me2'

data <- rbind(data_H3k9ac, data_H3k27ac, data_H3kcr, data_H3k4m2)
data$histone <- factor(data$histone, levels=c("H3K9ac", "H3K27ac", "Kcr", "H3K4me2"), ordered=T)
data$log2fc <- 0-data$log2fc

d1 <- data[data$histone=='H3K9ac',c("chr","start","end","type","interval","log2fc")]

interval_gene <- read.table('H3k9ac.consensus_peaks.annotatePeaks.txt', sep="\t", header=T)[,c(1,11)]
colnames(interval_gene) <- c("interval","Gene")
d1 <- merge(d1, interval_gene)

gene_log2fc <- read.table('LC_logfc.tsv', sep="\t", header=F)
colnames(gene_log2fc) <- c("Gene","TPM_log2fc")
gene_deg <- read.table('H0_H24_DEG.txt', sep="\t", header=T)
gene_info <- merge(gene_log2fc,gene_deg)

d1 <- merge(d1, gene_info)
str(d1)

d2 <- data[data$histone=='H3K4me2',c("chr","start","end","type","interval","log2fc")]

interval_gene <- read.table('H3k4m2.consensus_peaks.annotatePeaks.txt', sep="\t", header=T)[,c(1,11)]
colnames(interval_gene) <- c("interval","Gene")
d2 <- merge(d2, interval_gene)

gene_log2fc <- read.table('LC_logfc.tsv', sep="\t", header=F)
colnames(gene_log2fc) <- c("Gene","TPM_log2fc")
gene_deg <- read.table('H0_H24_DEG.txt', sep="\t", header=T)
gene_info <- merge(gene_log2fc,gene_deg)

d2 <- merge(d2, gene_info)
str(d2)

#d2 %>% count(type,DEG)

#pdf("NannoEpi_HiC_AB_HCLC_H3K4me2.pdf", useDingbats=FALSE, width=8, height=2.5, family="Arial", onefile=F)

pd1 <- ggscatter(d1, x="log2fc", y="TPM_log2fc", color="black", shape=21, size=0.3, facet.by="type",nrow=1,
          add="reg.line",  add.params=list(color = "blue", fill = "lightgray"),
          panel.labs = list(type = c("A2B", "A2A", "B2A", "B2B")),
          conf.int = TRUE, cor.coef = TRUE, 
          cor.coeff.args = list(method = "pearson", label.x = -3, label.sep = "\n")) +
    theme_bw() + theme(text=element_text(size=11, family = "Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("H3K9ac Peak log2FC") + ylab("Gene log2FC")

#dev.off()

#q()

pdf("NannoEpi_HiC_AB_HCLC_histone.pdf", useDingbats=FALSE, width=4, height=5, family="Arial", onefile=F)

ggviolin(data, x="type", y="log2fc", fill="type", palette="Paired", size=0.4, width=0.7) + 
     geom_boxplot(color="black", fill = "white", width=0.2, outlier.size=0.3) +
     facet_wrap(vars(histone), nrow=2, scales="free_y") + 
  scale_x_discrete(labels=c("TRUE_FALSE"="A2B","FALSE_TRUE"="B2A", "TRUE_TRUE"="A2A", "FALSE_FALSE"="B2B")) +
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Compartment changes") + ylab("Histone peaks (log2fc)") +
  stat_compare_means(comparisons=cmps2, method="wilcox.test",aes(labels=..p.signif..))

dev.off()

