#!/usr/bin/env Rscript

library(clusterProfiler)
library(enrichplot)

ann <- read.table("IMET1v2_genes.go.tsv",header=F,sep="\t",stringsAsFactors=F)

t2g <- data.frame(term=ann$V2,gene=ann$V1) 
gm <- buildGOmap(t2g)
tm <- go2term(gm$GO)
#tm2 <- go2term(t2g$term)

des <- read.table("DPNs_DEGs_Up.txt",header=F,stringsAsFactors=F)

res <- enricher(des$V1,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)
#res <- enricher(des$V1,TERM2GENE=t2g,TERM2NAME=tm2,pvalueCutoff = 0.05)
res

gl <- read.table("LC_logfc.tsv",header=F,sep="\t",stringsAsFactors=F)
geneList <- gl$V2
names(geneList) <- gl$V1
geneList = sort(geneList, decreasing = TRUE)

ego <- GSEA(geneList,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)

library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

pdf("GO_Up.pdf", useDingbats = FALSE, width = 12, height = 12, family="Arial")

a <- barplot(res,showCategory=20,x='GeneRatio') + theme_bw() +
    theme(text=element_text(size=11, family="Arial"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

b <- dotplot(res,showCategory=20) + theme_bw() +
    theme(text=element_text(size=11, family="Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

c <- cnetplot(res, foldChange=geneList, showCategory=20) + theme_bw() +
    theme(text=element_text(size=11, family="Arial"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

#d <- ridgeplot(ego)
edo <- pairwise_termsim(res)
d <- treeplot(edo)

ggarrange(a,b,c,d,nrow=2,ncol=2)

dev.off()

