#!/usr/bin/env Rscript

library(AnnotationForge)

IMET1v2_chromosome <- read.table("IMET1v2_chromosome.txt", sep="\t", header=F, stringsAsFactors=F)
colnames(IMET1v2_chromosome) <- c("GID", "CHROMOSOME")

IMET1v2_gene <- data.frame(GID=IMET1v2_chromosome$GID, SYMBOL=IMET1v2_chromosome$GID, GENENAME=IMET1v2_chromosome$GID)

IMET1v2_go <- read.table("IMET1v2_genes.go.tsv", sep="\t", header=F, stringsAsFactors=F)
colnames(IMET1v2_go) <- c("GID","GO")
IMET1v2_go$EVIDENCE <- "IEA"

makeOrgPackage(gene_info=IMET1v2_gene[!duplicated(IMET1v2_gene),], 
               chromosome=IMET1v2_chromosome[!duplicated(IMET1v2_chromosome),], 
               go=IMET1v2_go[!duplicated(IMET1v2_go),],
               version="0.1",
               maintainer="Yanhai Gong <gongyh@qibebt.ac.cn>",
               author="Yanhai Gong <gongyh@qibebt.ac.cn>",
               outputDir = ".",
               tax_id="1333499",
               genus="Nannochloropsis",
               species="oceanica",
               goTable="go")


install.packages("./org.Noceanica.eg.db", repos=NULL)
