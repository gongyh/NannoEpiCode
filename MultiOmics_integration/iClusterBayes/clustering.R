#!/usr/bin/env Rscript

# 加载 iClusterBayes
library(iClusterPlus)

rna_seq <- read.table("H0H24.TMM.TPM.matrix", header = TRUE, row.names = 1, sep = "\t")
h3k27ac <- read.table("H3k27ac/H3K27ac_normalized_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
h3k4me2 <- read.table("H3k4m2/H3K4me2_normalized_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
h3k9ac <- read.table("H3k9ac/H3K9ac_normalized_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
kcr <- read.table("H3kcr/Kcr_normalized_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
mnase_seq <- read.table("mnase_normalized_matrix.txt", header = TRUE, row.names = 1, sep = "\t")

samples <- c("H0_R1", "H0_R2", "H24_R1", "H24_R2")

# 数据矩阵
dt1 <- t(rna_seq)[samples,]   # RNA-Seq 数据
dt2 <- t(h3k4me2)[samples,]   # H3K4me2 ChIP-Seq 数据
dt3 <- t(h3k27ac)[samples,]   # H3K27ac ChIP-Seq 数据
dt4 <- t(kcr)[samples,]       # Kcr ChIP-Seq 数据
dt5 <- t(h3k9ac)[samples,]    # H3K9ac ChIP-Seq 数据
dt6 <- t(mnase_seq)[samples,] # MNase-Seq 数据

getHVFs <- function(dt) {
  # 计算CV
  gene_mean <- colMeans(dt)
  gene_sd <- apply(dt, 2, sd)
  gene_cv <- gene_sd / gene_mean

  # 按CV排序
  gene_cv_sorted <- sort(gene_cv, decreasing = TRUE)

  # 计算累计CV贡献
  cumulative_cv <- cumsum(gene_cv_sorted) / sum(gene_cv_sorted)

  # 选择解释 90% CV的特征
  n_features <- which(cumulative_cv >= 0.90)[1]
  hvg <- names(gene_cv_sorted)[1:n_features]

  # 返回新矩阵
  dt[, hvg]
}

dt1v <- getHVFs(dt1)
dt2v <- getHVFs(dt2*ncol(dt2)/10000)
dt3v <- getHVFs(dt3*ncol(dt3)/10000)
dt4v <- getHVFs(dt4*ncol(dt4)/10000)
dt5v <- getHVFs(dt5*ncol(dt5)/10000)
dt6v <- getHVFs(dt6*ncol(dt6)/10000)

# 指定每个数据类型的分布类型
#types <- c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian")
types <- c("gaussian", "poisson", "poisson", "poisson", "poisson", "poisson")

# 运行 iClusterBayes
set.seed(211)
fit <- iClusterBayes(dt1=dt1v, dt2=dt2v, dt3=dt3v, dt4=dt4v,
                     dt5=dt5v, dt6=dt6v, type=types, K=1)
save(fit, file="iClusterBayes.fit.Rdata")
#load("iClusterBayes.fit.Rdata")
str(fit)

# 查看聚类结果
print(fit$clusters)

sig_genes <- colnames(dt1v)[which(abs(fit$beta.pp[[1]])>0.5)]
write.table(sig_genes, file="iClusterBayes.sig.genes.csv", sep=",", quote=F, row.names=F, col.names=F)

h3k4me2_anno <- read.table("H3k4m2/H3k4m2.consensus_peaks.boolean.annotatePeaks.txt", header = TRUE, 
                           row.names = NULL, sep = "\t")
sig_h3k4me2 <- colnames(dt2v)[which(abs(fit$beta.pp[[2]])>0.5)]
sig_h3k4me2_g <- h3k4me2_anno[h3k4me2_anno$interval_id %in% sig_h3k4me2, c("interval_id","Entrez.ID")]
write.table(sig_h3k4me2_g, file="iClusterBayes.sig.h3k4me2.csv", sep=",", quote=F, row.names=F, col.names=F)
sig_h3k4me2_gene <- unique(sig_h3k4me2_g$Entrez.ID)

h3k27ac_anno <- read.table("H3k27ac/H3k27ac.consensus_peaks.boolean.annotatePeaks.txt", header = TRUE,
                           row.names = NULL, sep = "\t")
sig_h3k27ac <- colnames(dt3v)[which(abs(fit$beta.pp[[3]])>0.5)]
sig_h3k27ac_g <- h3k27ac_anno[h3k27ac_anno$interval_id %in% sig_h3k27ac, c("interval_id","Entrez.ID")]
write.table(sig_h3k27ac_g, file="iClusterBayes.sig.h3k27ac.csv", sep=",", quote=F, row.names=F, col.names=F)
sig_h3k27ac_gene <- unique(sig_h3k27ac_g$Entrez.ID)

kcr_anno <- read.table("H3kcr/H3kcr.consensus_peaks.boolean.annotatePeaks.txt", header = TRUE,
                       row.names = NULL, sep = "\t")
sig_kcr <- colnames(dt4v)[which(abs(fit$beta.pp[[4]])>0.5)]
sig_kcr_g <- kcr_anno[kcr_anno$interval_id %in% sig_kcr, c("interval_id","Entrez.ID")]
write.table(sig_kcr_g, file="iClusterBayes.sig.kcr.csv", sep=",", quote=F, row.names=F, col.names=F)
sig_kcr_gene <- unique(sig_kcr_g$Entrez.ID)

h3k9ac_anno <- read.table("H3k9ac/H3k9ac.consensus_peaks.boolean.annotatePeaks.txt", header = TRUE,
                          row.names = NULL, sep = "\t")
sig_h3k9ac <- colnames(dt5v)[which(abs(fit$beta.pp[[5]])>0.5)]
sig_h3k9ac_g <- h3k9ac_anno[h3k9ac_anno$interval_id %in% sig_h3k9ac, c("interval_id","Entrez.ID")]
write.table(sig_h3k9ac_g, file="iClusterBayes.sig.h3k9ac.csv", sep=",", quote=F, row.names=F, col.names=F)
sig_h3k9ac_gene <- unique(sig_h3k9ac_g$Entrez.ID)

mnase_anno <- read.table("mnase.consense_peaks.annotatePeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
colnames(mnase_anno)[1] <- "interval_id"
sig_mnase <- colnames(dt6v)[which(abs(fit$beta.pp[[6]])>0.5)]
sig_mnase_g <- mnase_anno[mnase_anno$interval_id %in% sig_mnase, c("interval_id","Entrez.ID")]
write.table(sig_mnase_g, file="iClusterBayes.sig.mnase.csv", sep=",", quote=F, row.names=F, col.names=F)
sig_mnase_gene <- unique(sig_mnase_g$Entrez.ID)

sig_epi_gene <- unique(c(sig_h3k4me2_gene,sig_h3k27ac_gene,sig_kcr_gene,sig_h3k9ac_gene,sig_mnase_gene))

## draw overlap with DEGs, DHMs, or DPNs
DEGs <- read.table("LC_HC_DE_info.txt", header = TRUE, row.names = NULL, sep = "\t")[,"gene"]
DHM_H3K4me2 <- read.table("H3k4m2/H3k4m2_DEPeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
DHM_H3K27ac <- read.table("H3k27ac/H3k27ac_DEPeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
DHM_H3K9ac <- read.table("H3k9ac/H3k9ac_DEPeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
DHM_Kcr <- read.table("H3kcr/H3kcr_DEPeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
DPNs <- read.table("DPNs_smt_consensus.bed", header = F, row.names = NULL, sep = "\t")
library(VennDiagram)
library(gridExtra)
library(ggpubr)
pdf("Sig_genes1.pdf", width=6, height=4)
colors <- get_palette("npg", 4)
p1 <- venn.diagram(list(DEGs=DEGs, Sig_genes=sig_genes), force.unique=T, fill = c(colors[1],colors[2]), 
                   filename=NULL, lwd=1, scaled=T, alpha=0.6, cat.col=c(colors[1],colors[2]), margin=0.1)
p2 <- venn.diagram(list(DHM_H3K4me2=DHM_H3K4me2[DHM_H3K4me2$DEPeak != "Not", "Peak"], Sig_H3K4me2=sig_h3k4me2_g$interval_id), 
                   force.unique=T, filename=NULL, fill = c(colors[1],colors[2]), lwd=1, margin=0.1,
                   scaled = TRUE, alpha=0.6, cat.col = c(colors[1],colors[2]))
p3 <- venn.diagram(list(DHM_H3K9ac=DHM_H3K9ac[DHM_H3K9ac$DEPeak != "Not", "Peak"], Sig_H3K9ac=sig_h3k9ac_g$interval_id), 
                   filename=NULL, force.unique=T, fill = c(colors[1],colors[2]), lwd=1, margin=0.1,
                   scaled = TRUE, alpha=0.6, cat.col = c(colors[1],colors[2]))
p4 <- venn.diagram(list(DHM_H3K27ac=DHM_H3K27ac[DHM_H3K27ac$DEPeak != "Not", "Peak"], Sig_H3K27ac=sig_h3k27ac_g$interval_id), 
                   filename=NULL, force.unique=T, fill = c(colors[1],colors[2]), lwd=1, margin=0.1,
                   scaled = TRUE, alpha=0.6, cat.col = c(colors[1],colors[2]))
p5 <- venn.diagram(list(DHM_Kcr=DHM_Kcr[DHM_Kcr$DEPeak != "Not", "Peak"], Sig_Kcr=sig_kcr_g$interval_id), filename=NULL, 
                   force.unique=T, fill = c(colors[1],colors[2]), lwd=1, margin=0.1,
                   scaled = TRUE, alpha=0.6, cat.col = c(colors[1],colors[2]))
p6 <- venn.diagram(list(DPN_smt=paste0(DPNs[,1],":",DPNs[,2],"-",DPNs[,3]), Sig_DPN=sig_mnase_g$interval_id), filename=NULL, 
                   force.unique=T, fill = c(colors[1],colors[2]), lwd=1, margin=0.1,
                   scaled = TRUE, alpha=0.6, cat.col = c(colors[1],colors[2]))
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)
dev.off()

x <- list(H3K4me2=sig_h3k4me2_gene, H3K9ac=sig_h3k9ac_gene, H3K27ac=sig_h3k27ac_gene, Kcr=sig_kcr_gene)
pdf("Sig_genes3.pdf", width=4, height=4)
p7 <- venn.diagram(x, filename=NULL, force.unique=T, fill=colors, lwd=1, margin=0.1, scaled=T, alpha=0.6, cat.col=colors)
grid.draw(p7)
dev.off()

y <- list(mRNA=sig_genes, Nucleosome=sig_mnase_gene)
pdf("Sig_genes2.pdf", width=2, height=2)
p8 <- venn.diagram(y, filename=NULL, force.unique=T, fill = c(colors[1],colors[2]), lwd=1, margin=0.1,
             scaled = TRUE, alpha=0.6, cat.col = c(colors[1],colors[2]))
grid.draw(p8)
dev.off()

