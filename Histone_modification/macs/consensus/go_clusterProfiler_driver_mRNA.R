#!/usr/bin/env Rscript

setwd('./')

library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(DOSE)
library(dplyr)
library(org.Noceanica.eg.db)
library(extrafont)
loadfonts()
library(ggpubr)

print("------ sig genes -----")
dep_9ac <- read.table("iClusterBayes.sig.genes.csv",header=F,sep=",",stringsAsFactors=F)
dep_9ac_genes <- dep_9ac$V1

dep_27ac <- read.table("H0_H24_DEG.txt",header=T,sep="\t",stringsAsFactors=F)
logfc <- read.table("VLC24vsH0_logFC.txt",header=T,sep="\t",stringsAsFactors=F)
dep_27ac_genes <- dep_27ac[, "Gene"] #all genes
#dep_27ac_genes2 <- dep_27ac[dep_27ac$DEG!="Not", "Gene"] #DEGs
#dep_27ac_genes <- setdiff(dep_9ac_genes,  dep_27ac_genes2) #sig nonDEGs
gene_log2fc <- logfc[logfc$Gene %in% dep_27ac_genes,]
gene_list <- gene_log2fc$logFC
names(gene_list) <- gene_log2fc$Gene
gene_list_sorted <- sort(gene_list, decreasing=T)

# 进行 GSEA 分析
gsea_result1 <- gseGO(
  geneList = gene_list_sorted,  # 排序的基因列表
  ont = "BP",            # GO 子类别：BP（生物过程）、MF（分子功能）、CC（细胞组分）
  OrgDb = org.Noceanica.eg.db,  # 物种注释数据库
  keyType = "GID",    # 基因 ID 类型（如 SYMBOL、ENTREZID 等）
  exponent = 1,          # 权重指数（默认为 1）
  minGSSize = 1,        # 最小基因集大小
  maxGSSize = 1000,       # 最大基因集大小
  pvalueCutoff = 0.05,   # p 值阈值
  pAdjustMethod = "BH",  # p 值校正方法
  nPermSimple = 10000,
  eps = 0,
  verbose = FALSE        # 是否显示详细信息
)

gsea_BP <- simplify(gsea_result1)

gsea_result2 <- gseGO(
  geneList = gene_list_sorted,  # 排序的基因列表
  ont = "MF",            # GO 子类别：BP（生物过程）、MF（分子功能）、CC（细胞组分）
  OrgDb = org.Noceanica.eg.db,  # 物种注释数据库
  keyType = "GID",    # 基因 ID 类型（如 SYMBOL、ENTREZID 等）
  exponent = 1,          # 权重指数（默认为 1）
  minGSSize = 1,        # 最小基因集大小
  maxGSSize = 1000,       # 最大基因集大小
  pvalueCutoff = 0.05,   # p 值阈值
  pAdjustMethod = "BH",  # p 值校正方法
  nPermSimple = 10000,
  eps = 0,
  verbose = FALSE        # 是否显示详细信息
)

gsea_MF <- simplify(gsea_result2)

gsea_result3 <- gseGO(
  geneList = gene_list_sorted,  # 排序的基因列表
  ont = "CC",            # GO 子类别：BP（生物过程）、MF（分子功能）、CC（细胞组分）
  OrgDb = org.Noceanica.eg.db,  # 物种注释数据库
  keyType = "GID",    # 基因 ID 类型（如 SYMBOL、ENTREZID 等）
  exponent = 1,          # 权重指数（默认为 1）
  minGSSize = 1,        # 最小基因集大小
  maxGSSize = 1000,       # 最大基因集大小
  pvalueCutoff = 0.05,   # p 值阈值
  pAdjustMethod = "BH",  # p 值校正方法
  nPermSimple = 10000,
  eps = 0,
  verbose = FALSE        # 是否显示详细信息
)

gsea_CC <- simplify(gsea_result3)

# 查看结果

data <- rbind(gsea_BP@result, gsea_MF@result, gsea_CC@result)
data$Ontology <- Ontology(GOTERM[data$ID])

write.csv(data, file="GO_GSEA.csv")

core_sig <- apply(data, 1, function(row) {
  core_genes <- strsplit(row[11], split = "/")[[1]]
  core_genes_DEG <- setdiff(core_genes, dep_27ac[dep_27ac$DEG !="Not","Gene"])
  core_genes_Sig <- intersect(core_genes_DEG, dep_9ac_genes)
  core_genes_Sig
})

write.table(core_sig[["GO:0005622"]], file="GO0005622_Sig.txt", quote=F, row.names=F, col.names=F)

## plot three panels with the same height

library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

#pdf("GO_GSEA4.pdf", useDingbats = FALSE, width = 8, height = 2.5, family="Arial", onefile=F)

## 1. separate genes into 4 groups (Driver_DEG, Driver_nDEG, nDriver_DEG, nDriver_nDEG), 
##    then draw the percentage of epigenetic regulated
all_gene <- dep_27ac_genes
sig_gene <- dep_9ac_genes
nsig_gene <- setdiff(all_gene, sig_gene)
de_gene <- dep_27ac[dep_27ac$DEG!="Not", "Gene"]
sig_de <- intersect(sig_gene, de_gene)
sig_nde <- setdiff(sig_gene, de_gene)
nsig_de <- intersect(nsig_gene, de_gene)
nsig_nde <- setdiff(nsig_gene, de_gene)

# get all epi-regulated genes
DHM_H3K4me2 <- read.table("H3k4m2/H3k4m2_DEPeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
DHMG_H3K4me2 <- DHM_H3K4me2[DHM_H3K4me2$DEPeak != "Not", "Gene"]
DHM_H3K27ac <- read.table("H3k27ac/H3k27ac_DEPeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
DHMG_H3K27ac <- DHM_H3K27ac[DHM_H3K27ac$DEPeak != "Not", "Gene"]
DHM_H3K9ac <- read.table("H3k9ac/H3k9ac_DEPeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
DHMG_H3K9ac <- DHM_H3K9ac[DHM_H3K9ac$DEPeak != "Not", "Gene"]
DHM_Kcr <- read.table("H3kcr/H3kcr_DEPeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
DHMG_Kcr <- DHM_Kcr[DHM_Kcr$DEPeak != "Not", "Gene"]
DPNs <- read.table("DPNs_smt_consensus.bed", header = F, row.names = NULL, sep = "\t")
DPN_smt <- paste0(DPNs[,1],":",DPNs[,2],"-",DPNs[,3])
mnase_anno <- read.table("mnase.consense_peaks.annotatePeaks.txt", header = TRUE, row.names = NULL, sep = "\t")
colnames(mnase_anno)[1] <- "interval_id"
DPNG_smt <- mnase_anno[mnase_anno$interval_id %in% DPN_smt, "Entrez.ID"]
epi_gene <- unique(c(DHMG_H3K4me2,DHMG_H3K27ac,DHMG_H3K9ac,DHMG_Kcr,DPNG_smt))

# get epi percent
sig_de_epi <- intersect(sig_de, epi_gene)
sig_nde_epi <- intersect(sig_nde, epi_gene)
nsig_de_epi <- intersect(nsig_de, epi_gene)
nsig_nde_epi <- intersect(nsig_nde, epi_gene)
sig_de_epi_ratio <- 100 * length(sig_de_epi) / length(sig_de)
sig_nde_epi_ratio <- 100 * length(sig_nde_epi) / length(sig_nde)
nsig_de_epi_ratio <- 100 * length(nsig_de_epi) / length(nsig_de)
nsig_nde_epi_ratio <- 100 * length(nsig_nde_epi) / length(nsig_nde)

ratio_df <- data.frame(type=c("Sig_DEG","Sig_nDEG","nSig_DEG","nSig_nDEG"), 
  ratio=c(sig_de_epi_ratio,sig_nde_epi_ratio,nsig_de_epi_ratio,nsig_nde_epi_ratio))

pdf("GO_GSEA4.pdf", useDingbats = FALSE, width = 8, height = 2.5, family="Arial", onefile=F)

p1 <- ggbarplot(ratio_df, x="type", y="ratio", fill="type", color="type", palette="npg", 
                label=TRUE, label.pos="out", lab.size=2, lab.nb.digits=1, width=0.5) +
    theme_bw() + xlab(NULL) + ylab("Epigenetically regulated (%)") +
    theme(text=element_text(size=9, family="Arial"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", angle = 30, hjust = 1),
        axis.text.y = element_text(color="black"))

## 2. plot GSEA results, classify core_enrichment genes into the four groups.

data$Description <- stringr::str_to_sentence(data$Description)

data$ce4 <- sapply(data$core_enrichment, function(ce) {
  core_genes <- strsplit(ce, split = "/")[[1]]
  cg_sig_de <- intersect(core_genes, sig_de)
  cg_sig_nde <- intersect(core_genes, sig_nde)
  cg_nsig_de <- intersect(core_genes, nsig_de)
  cg_nsig_nde <- intersect(core_genes, nsig_nde)
  paste(length(cg_sig_de),length(cg_sig_nde),
        length(cg_nsig_de),length(cg_nsig_nde),sep=",")
})

col3 <- get_palette(palette = "default", 3)
data <- data %>%
  mutate(color_Ontology = case_when(
    Ontology == "MF" ~ col3[3],
    Ontology == "CC" ~ col3[2],
    Ontology == "BP" ~ col3[1]
  ))

p2 <- ggbarplot(data, x="ID", y="qvalue", fill="NES", color="NES", group="Ontology", sort.val="none",
                sort.by.groups=T, width=0.1, orientation = "horiz") + 
          scale_x_discrete(expand=expansion(add=c(0.5,1))) +
          scale_y_continuous(trans=reverselog_trans(10), limits=c(1,1e-37), expand = c(0, 0.1),
                             breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=9e-1,label=Description), hjust = 0, vjust=-0.65, family="Arial", size=2)+
    geom_text(aes(y=1e-30,label=ce4), hjust = 0, vjust=0.5, family="Arial", size=2)+
    scale_fill_viridis_c(option = "rocket") + scale_color_viridis_c(option = "rocket") +
    theme_bw() + xlab(NULL) + ylab(NULL) +
    theme(text=element_text(size=9, family="Arial"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"), legend.key.height = unit(0.2, "cm"),
        legend.box.spacing = unit(0.05, "cm"), legend.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color=data$color_Ontology),
        axis.text.x = element_text(color="black"))

## 3. visualize Driver_nDEG for core_sig[["GO:0005622"]]

cc3 <- gsearank(gsea_CC, 3, title = "Intracellular anatomical structure", output = "data")
cc3_logfc <- merge(cc3, logfc, by=1)
cc3_core <- cc3_logfc[cc3_logfc[,"core enrichment"]=="YES",]
cc3_core$type <- sapply(cc3_core$gene, function(g) {
  if (g %in% sig_de) {
    "Sig_DEG"
  } else if(g %in% sig_nde) {
    "Sig_nDEG"
  } else if(g %in% nsig_de) {
    "nSig_DEG"
  } else if(g %in% nsig_nde) {
    "nSig_nDEG"
  }
})

cc3_core$type <- factor(cc3_core$type, levels=c("Sig_DEG","Sig_nDEG","nSig_DEG","nSig_nDEG"), ordered=T)
p3 <- ggscatter(cc3_core, x="logFC", y="running ES", color="type",shape="type",
                palette="npg",rug=F,size=0.5) + labs(title="Intracellular anatomical structure") +
  annotate("text", x=-1.5, y=-0.05, label="NES: -1.54\nQ: 1.04e-11", size=2.5) +
  theme_bw() + xlab("Gene expression (log2FC)") + ylab("Running ES") +
    theme(text=element_text(size=9, family="Arial"),
        plot.title = element_text(color="black", size=9),
        legend.position = c(0.25,0.35), legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black"))

ggarrange(p1,p2,p3, widths=c(1,2,1.5), align="h", nrow=1)

dev.off()

