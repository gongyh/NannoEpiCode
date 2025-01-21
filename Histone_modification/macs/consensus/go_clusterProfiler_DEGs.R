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

print("------ DEGs -----")
dep_27ac <- read.table("H0_H24_DEG.txt",header=T,stringsAsFactors=F)
dep_27ac_genes <- dep_27ac[dep_27ac$DEG != "Not", "Gene"]
res <- enrichGO(dep_27ac_genes, 'org.Noceanica.eg.db', ont="ALL", pvalueCutoff=0.01, keyType = "GID")
#res

res_bp <- enrichGO(dep_27ac_genes, 'org.Noceanica.eg.db', ont="BP", pvalueCutoff=0.01, keyType = "GID")
res_bp2 <- simplify(res_bp)
res_mf <- enrichGO(dep_27ac_genes, 'org.Noceanica.eg.db', ont="MF", pvalueCutoff=0.01, keyType = "GID")
res_mf2 <- simplify(res_mf)
res_cc <- enrichGO(dep_27ac_genes, 'org.Noceanica.eg.db', ont="CC", pvalueCutoff=0.01, keyType = "GID")
res_cc2 <- simplify(res_cc)

data2 <- rbind(res_bp2@result, res_mf2@result, res_cc2@result)
data2$Ontology <- Ontology(GOTERM[data2$ID])

print("------ Done! -----")

pdf("GO_DEGs_bar.pdf", useDingbats = FALSE, width = 3, height = 3, family="Arial")

library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

data2$Description <- stringr::str_to_sentence(data2$Description)

str(data2)

ggbarplot(data2, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
          sort.by.groups=T, width=0.1, orientation = "horiz") + scale_x_discrete(expand = expansion(add = c(0.5,1))) +
          scale_y_continuous(trans=reverselog_trans(10),limits=c(1,1e-12), expand = c(0, 0.1),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=1,label=Description), hjust = 0, vjust=-1, family="Arial", size=2)+
    geom_text(aes(y=1e-11,label=Count), hjust = 0, vjust=0.5, family="Arial", size=2)+
    guides(fill = guide_legend(reverse = TRUE), colour = "none") +
    theme_bw() + xlab("") + ylab("Q value") + 
    theme(text=element_text(size=9, family="Arial"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

dev.off()

