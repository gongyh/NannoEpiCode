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

print("------ H3K9ac -----")
dep_9ac <- read.table("H3k9ac/H3k9ac_DEPeaks.txt",header=T,stringsAsFactors=F)
dep_9ac_genes <- dep_9ac[dep_9ac$DEPeak != "Not", "Gene"]

res_bp <- enrichGO(dep_9ac_genes, 'org.Noceanica.eg.db', ont="BP", pvalueCutoff=0.01, keyType = "GID")
res_bp2 <- simplify(res_bp)
res_mf <- enrichGO(dep_9ac_genes, 'org.Noceanica.eg.db', ont="MF", pvalueCutoff=0.01, keyType = "GID")
res_mf2 <- simplify(res_mf)
res_cc <- enrichGO(dep_9ac_genes, 'org.Noceanica.eg.db', ont="CC", pvalueCutoff=0.01, keyType = "GID")
res_cc2 <- simplify(res_cc)

#data1 <- rbind(res_bp2@result, res_mf2@result, res_cc2@result)
#data1$Ontology <- Ontology(GOTERM[data1$ID])
#data1$histone <- 'H3K9ac'

res <- enrichGO(dep_9ac_genes, 'org.Noceanica.eg.db', ont="ALL", pvalueCutoff=0.01, keyType = "GID")
#res

print("------ H3K27ac -----")
dep_27ac <- read.table("H3k27ac/H3k27ac_DEPeaks.txt",header=T,stringsAsFactors=F)
dep_27ac_genes <- dep_27ac[dep_27ac$DEPeak != "Not", "Gene"]
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
data2$histone <- 'H3K27ac'

print("------ Kcr -----")
dep_cr <- read.table("H3kcr/H3kcr_DEPeaks.txt",header=T,stringsAsFactors=F)
dep_cr_genes <- dep_cr[dep_cr$DEPeak != "Not", "Gene"]
res <- enrichGO(dep_cr_genes, 'org.Noceanica.eg.db', ont="ALL", pvalueCutoff=0.01, keyType = "GID")
#res

res_bp <- enrichGO(dep_cr_genes, 'org.Noceanica.eg.db', ont="BP", pvalueCutoff=0.01, keyType = "GID")
res_bp2 <- simplify(res_bp)
res_mf <- enrichGO(dep_cr_genes, 'org.Noceanica.eg.db', ont="MF", pvalueCutoff=0.01, keyType = "GID")
res_mf2 <- simplify(res_mf)
res_cc <- enrichGO(dep_cr_genes, 'org.Noceanica.eg.db', ont="CC", pvalueCutoff=0.01, keyType = "GID")
res_cc2 <- simplify(res_cc)

data3 <- rbind(res_bp2@result, res_mf2@result, res_cc2@result)
data3$Ontology <- Ontology(GOTERM[data3$ID])
data3$histone <- 'Kcr'

print("------ H3K4me2 -----")
dep_4m2 <- read.table("H3k4m2/H3k4m2_DEPeaks.txt",header=T,stringsAsFactors=F)
dep_4m2_genes <- dep_4m2[dep_4m2$DEPeak != "Not", "Gene"]
res <- enrichGO(dep_4m2_genes, 'org.Noceanica.eg.db', ont="ALL", pvalueCutoff=0.01, keyType = "GID")
#res
res_bp <- enrichGO(dep_4m2_genes, 'org.Noceanica.eg.db', ont="BP", pvalueCutoff=0.01, keyType = "GID")
res_bp2 <- simplify(res_bp)
res_mf <- enrichGO(dep_4m2_genes, 'org.Noceanica.eg.db', ont="MF", pvalueCutoff=0.01, keyType = "GID")
res_mf2 <- simplify(res_mf)
res_cc <- enrichGO(dep_4m2_genes, 'org.Noceanica.eg.db', ont="CC", pvalueCutoff=0.01, keyType = "GID")
res_cc2 <- simplify(res_cc)

data4 <- rbind(res_bp2@result, res_mf2@result, res_cc2@result)
data4$Ontology <- Ontology(GOTERM[data4$ID])
data4$histone <- 'H3K4me2'

print("------ Done! -----")

library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

data2$Description <- stringr::str_to_sentence(data2$Description)
data3$Description <- stringr::str_to_sentence(data3$Description)
data4$Description <- stringr::str_to_sentence(data4$Description)

write.csv(data2, file="GO_DHM_H3K27ac.csv")
write.csv(data3, file="GO_DHM_Kcr.csv")
write.csv(data4, file="GO_DHM_H3K4me2.csv")

q()

pdf("GO_DHM_bar.pdf", useDingbats = FALSE, width = 10, height = 2.5, family="Arial")

p2 <- ggbarplot(data2, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
          sort.by.groups=T, width=0.1, orientation = "horiz") + scale_x_discrete(expand = expansion(add = c(0.5,1))) +
          scale_y_continuous(trans=reverselog_trans(10),limits=c(1,1e-10), expand = c(0, 0.1),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=1,label=Description), hjust = 0, vjust=-1, family="Arial", size=2)+
    geom_text(aes(y=1e-9,label=Count), hjust = -0.5, vjust=0.5, family="Arial", size=2)+
    guides(fill = guide_legend(reverse = TRUE), colour = "none") +
    theme_bw() + xlab("") + ylab("Q value") + 
    theme(text=element_text(size=9, family="Arial"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

p3 <- ggbarplot(data3, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
          sort.by.groups=T, width=0.1, orientation = "horiz") + scale_x_discrete(expand = expansion(add = c(0.5,1))) +
          scale_y_continuous(trans=reverselog_trans(10),limits=c(1,1e-5), expand = c(0, 0.1),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=1,label=Description), hjust = 0, vjust=-1, family="Arial", size = 2)+
    geom_text(aes(y=3e-5,label=Count), hjust = -0.5, vjust=0.5, family="Arial", size = 2)+
    guides(fill = guide_legend(reverse = TRUE), colour = "none") +
    theme_bw() + xlab("") + ylab("Q value") +
    theme(text=element_text(size=9, family="Arial"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

p4 <- ggbarplot(data4, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
          sort.by.groups=T, width=0.1, orientation = "horiz") + scale_x_discrete(expand = expansion(add = c(0.5,1))) +
          scale_y_continuous(trans=reverselog_trans(10),limits=c(1,1e-20), expand = c(0, 0.1),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=1,label=Description), hjust = 0, vjust=-1, family="Arial", size = 2)+
    geom_text(aes(y=1e-18,label=Count), hjust = -0.5, vjust=0.5, family="Arial", size = 2)+
    guides(fill = guide_legend(reverse = TRUE), colour = "none") +
    theme_bw() + xlab("") + ylab("Q value") +
    theme(text=element_text(size=9, family="Arial"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

ggarrange(p2,p3,p4, nrow=1, common.legend=T, legend="top")

dev.off()

