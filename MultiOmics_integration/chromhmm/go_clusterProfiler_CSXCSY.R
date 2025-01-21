#!/usr/bin/env Rscript

library(clusterProfiler)
library(enrichplot)
library(GO.db)

ann <- read.table("IMET1v2_genes.go.tsv",header=F,sep="\t",stringsAsFactors=F)

t2g <- data.frame(term=ann$V2,gene=ann$V1) 
gm <- buildGOmap(t2g)
tm <- go2term(gm$GO)
#tm2 <- go2term(t2g$term)

des <- read.table("CS1_CS3_genes.txt",header=F,stringsAsFactors=F)
res13 <- enricher(des$V1,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)
res13_table <- res13@result
data13 <- res13_table[res13_table$qvalue<0.01,]
str(data13)

des <- read.table("CS1_CS4_genes.txt",header=F,stringsAsFactors=F)
res14 <- enricher(des$V1,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)
res14_table <- res14@result
data14 <- res14_table[res14_table$qvalue<0.01,]
str(data14)

des <- read.table("CS2_CS3_genes.txt",header=F,stringsAsFactors=F)
res23 <- enricher(des$V1,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)
res23_table <- res23@result
data23 <- res23_table[res23_table$qvalue<0.01,]
str(data23)

des <- read.table("CS3_CS1_genes.txt",header=F,stringsAsFactors=F)
res31 <- enricher(des$V1,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)
res31_table <- res31@result
data31 <- res31_table[res31_table$qvalue<0.01,]
str(data31)

des <- read.table("CS3_CS2_genes.txt",header=F,stringsAsFactors=F)
res32 <- enricher(des$V1,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)
res32_table <- res32@result
data32 <- res32_table[res32_table$qvalue<0.01,]
str(data32)

des <- read.table("CS4_CS1_genes.txt",header=F,stringsAsFactors=F)
res41 <- enricher(des$V1,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)
res41_table <- res41@result
data41 <- res41_table[res41_table$qvalue<0.01,]
str(data41)

des <- read.table("CS4_CS2_genes.txt",header=F,stringsAsFactors=F)
res42 <- enricher(des$V1,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)
res42_table <- res42@result
data42 <- res42_table[res42_table$qvalue<0.01,]
str(data42)

des <- read.table("CS5_CS2_genes.txt",header=F,stringsAsFactors=F)
res52 <- enricher(des$V1,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)
res52_table <- res52@result
data52 <- res52_table[res52_table$qvalue<0.01,]
str(data52)

q()

gl <- read.table("LC_logfc.tsv",header=F,sep="\t",stringsAsFactors=F)
geneList <- gl$V2
names(geneList) <- gl$V1
geneList = sort(geneList, decreasing = TRUE)

ego <- GSEA(geneList,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)

library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

pdf("GO_CSXCSY.pdf", useDingbats = FALSE, width = 5, height = 7.5, family="Arial")

res_table <- res@result
data <- res_table[res_table$qvalue<0.01,]
data$Ontology <- Ontology(GOTERM[data$ID])

ggbarplot(data, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
          sort.by.groups=T, width=0.1, orientation = "horiz") +
          scale_y_continuous(trans=reverselog_trans(10),limits=c(1,1e-15), expand = c(0, 0.1),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=1,label=Description), hjust = 0, vjust=-1, family="Arial", size = 3)+
    geom_text(aes(y=1e-14,label=Count), hjust = -0.5, vjust=0.5, family="Arial", size = 3)+
    theme_bw() + xlab("") + ylab("Q value") + 
    theme(text=element_text(size=11, family="Arial"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

dev.off()

