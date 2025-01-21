#!/usr/bin/env Rscript

setwd('./')

library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(DOSE)
library(dplyr)
library(org.Noceanica.eg.db)

des <- read.table("../DPNs_DEGs_Down.txt",header=F,stringsAsFactors=F)

res_bp <- enrichGO(des$V1, 'org.Noceanica.eg.db', ont="BP", pvalueCutoff=0.01, keyType = "GID")
res_bp2 <- simplify(res_bp)
res_mf <- enrichGO(des$V1, 'org.Noceanica.eg.db', ont="MF", pvalueCutoff=0.01, keyType = "GID")
res_mf2 <- simplify(res_mf)
res_cc <- enrichGO(des$V1, 'org.Noceanica.eg.db', ont="CC", pvalueCutoff=0.01, keyType = "GID")
res_cc2 <- simplify(res_cc)

res <- enrichGO(des$V1, 'org.Noceanica.eg.db', ont="ALL", pvalueCutoff=0.01, keyType = "GID")

gl <- read.table("../LC_logfc.tsv",header=F,sep="\t",stringsAsFactors=F)
geneList <- gl$V2
names(geneList) <- gl$V1
geneList = sort(geneList, decreasing = TRUE)

#ego <- GSEA(geneList,TERM2GENE=gm,TERM2NAME=tm,pvalueCutoff = 0.05)

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

list2df <- function(inputList) {
    # ldf <- lapply(1:length(inputList), function(i) {
    ldf <- lapply(seq_len(length(inputList)), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]])
    })

    do.call('rbind', ldf)
}

pdf("GO_DPNd_bar.pdf", useDingbats = FALSE, width = 5, height = 3, family="Arial")

res_table <- rbind(res_bp2@result, res_mf2@result, res_cc2@result)
data <- res_table[res_table$p.adjust<0.01,]
data$Ontology <- Ontology(GOTERM[data$ID])

#str(data)

# manually remove duplicates
dupGOs <- c("GO:0043604", "GO:0006518", "GO:0043043", "GO:0034645", "GO:0010467")
data <- data[! data$ID %in% dupGOs, ]

ggbarplot(data, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
          sort.by.groups=T, width=0.1, orientation = "horiz") +
          scale_y_continuous(trans=reverselog_trans(10),limits=c(1,1e-15), expand = c(0, 0.1),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=1,label=Description), hjust = 0, vjust=-1, family="Arial", size = 3)+
    geom_text(aes(y=1e-14,label=Count), hjust = -0.5, vjust=0.5, family="Arial", size = 3)+
    guides(fill = guide_legend(reverse = TRUE), colour = "none") +
    theme_bw() + xlab("") + ylab("Q value") + 
    theme(text=element_text(size=11, family="Arial"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

dev.off()

pdf("GO_DPNd_upset.pdf", useDingbats = FALSE, width = 5, height = 4, family="Arial")

df <- data
id <- df$ID
des <- df$Description
glist <- apply(data,1,function(x){strsplit(x["geneID"],"/")[[1]]},simplify=F)
#names(glist) <- des
d <- list2df(glist)
res <- tibble::tibble(Description = split(d[,1], d[,2]))
ggplot(res, aes_(x = ~Description)) + geom_bar() +
    geom_text(stat='count', aes(label=after_stat(count)), vjust=1, colour="white") +
    theme_dose(font.size = 12) +
    xlab(NULL) + ylab(NULL) +
    ggupset::scale_x_upset(order_by = "freq")

dev.off()

