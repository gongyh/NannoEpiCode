#!/usr/bin/env Rscript

library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(DOSE)

ann <- read.table("IMET1v2_genes.go.tsv",header=F,sep="\t",stringsAsFactors=F)

t2g <- data.frame(term=ann$V2,gene=ann$V1) 
gm <- NULL
if (file.exists('gm.rds')) {
  gm <- readRDS("gm.rds")
} else {
  gm <- buildGOmap(t2g)
  saveRDS(gm, "gm.rds")
}
tm <- go2term(gm$GO)
#tm2 <- go2term(t2g$term)

des <- read.table("DPNs_DEGs_Down.txt",header=F,stringsAsFactors=F)

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

list2df <- function(inputList) {
    # ldf <- lapply(1:length(inputList), function(i) {
    ldf <- lapply(seq_len(length(inputList)), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]])
    })

    do.call('rbind', ldf)
}

library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

pdf("GO_Down_bar.pdf", useDingbats = FALSE, width = 5, height = 7.5, family="Arial")

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

pdf("GO_Down_upset.pdf", useDingbats = FALSE, width = 10, height = 10, family="Arial")

edox <- pairwise_termsim(res)

#cnetplot(res, foldChange=geneList, showCategory=20, circular=T, colorEdge=T) + theme_bw() +
#    theme(text=element_text(size=11, family="Arial"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        axis.text.x = element_text(color="black"),
#        axis.text.y = element_text(color="black"))

#cnetplot(res, foldChange=geneList, showCategory=20, node_label="gene") + theme_bw() +
#    theme(text=element_text(size=11, family="Arial"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        axis.text.x = element_text(color="black"),
#        axis.text.y = element_text(color="black"))

#treeplot(edox)
#emapplot(edox)

#upsetplot(edox, n=20)
    n <- 20
    df <- as.data.frame(edox)
    id <- df$ID[1:n]
    des <- df$Description[1:n]
    glist <- geneInCategory(edox)[id]
    names(glist) <- des
    d <- list2df(glist)
    res <- tibble::tibble(Description = split(d[,1], d[,2]))
    ggplot(res, aes_(x = ~Description)) + geom_bar() +
        geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
        theme_dose(font.size = 12) +
        xlab(NULL) + ylab(NULL) +
        ggupset::scale_x_upset(order_by = "freq")


dev.off()

