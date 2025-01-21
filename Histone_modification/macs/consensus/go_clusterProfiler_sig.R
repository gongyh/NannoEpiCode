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
dep_9ac <- read.table("iClusterBayes.sig.h3k9ac.csv",header=F,sep=",",stringsAsFactors=F)
dep_9ac_genes <- dep_9ac$V2

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
dep_27ac <- read.table("iClusterBayes.sig.h3k27ac.csv",header=F,sep=",",stringsAsFactors=F)
dep_27ac_genes <- dep_27ac$V2
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
dep_cr <- read.table("iClusterBayes.sig.kcr.csv",header=F,sep=",",stringsAsFactors=F)
dep_cr_genes <- dep_cr$V2
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
dep_4m2 <- read.table("iClusterBayes.sig.h3k4me2.csv",header=F,sep=",",stringsAsFactors=F)
dep_4m2_genes <- dep_4m2$V2
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

write.csv(data2, file="GO_sig_H3K27ac.csv")
write.csv(data3, file="GO_sig_Kcr.csv")
write.csv(data4, file="GO_sig_H3K4me2.csv")

pdf("GO_Sig_bar.pdf", useDingbats = FALSE, width = 6, height = 2.5, family="Arial", onefile=F)

p2 <- ggbarplot(data2, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
          sort.by.groups=T, width=0.1, orientation = "horiz") + scale_x_discrete(expand = expansion(add = c(0.5,1))) +
          scale_y_continuous(trans=reverselog_trans(10), limits=c(1,1e-14), expand = c(0, 0.1),
                             breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=1,label=Description), hjust = 0, vjust=-0.65, family="Arial", size=2)+
    geom_text(aes(y=1e-13,label=Count), hjust = -0.5, vjust=0.5, family="Arial", size=2)+
    guides(fill = guide_legend(reverse = TRUE), colour = "none") +
    theme_bw() + xlab("") + ylab("Q value") + 
    theme(text=element_text(size=9, family="Arial"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

p3 <- ggbarplot(data3, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
          sort.by.groups=T, width=0.1, orientation = "horiz") + scale_x_discrete(expand = expansion(add = c(0.5,1))) +
          scale_y_continuous(trans=reverselog_trans(10),limits=c(1,1e-4), expand = c(0, 0.1),
                             breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=1,label=Description), hjust = 0, vjust=-0.65, family="Arial", size = 2)+
    geom_text(aes(y=2e-4,label=Count), hjust = -0.5, vjust=0.5, family="Arial", size = 2)+
    guides(fill = guide_legend(reverse = TRUE), colour = "none") +
    theme_bw() + labs(x=NULL,y=NULL) +
    theme(text=element_text(size=9, family="Arial"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

p4 <- ggbarplot(data4, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
          sort.by.groups=T, width=0.1, orientation = "horiz") + scale_x_discrete(expand = expansion(add = c(0.5,1))) +
          scale_y_continuous(trans=reverselog_trans(10),limits=c(1,1e-4), expand = c(0, 0.1),
                             breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x))) +
    geom_text(aes(y=1,label=Description), hjust = 0, vjust=-0.65, family="Arial", size = 2)+
    geom_text(aes(y=2e-4,label=Count), hjust = -0.5, vjust=0.5, family="Arial", size = 2)+
    guides(fill = guide_legend(reverse = TRUE), colour = "none") +
    theme_bw() + xlab("") + ylab("Q value") +
    theme(text=element_text(size=9, family="Arial"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

ggarrange(p2,ggarrange(p3,p4,ncol=1,heights=c(2,1.25),common.legend=T, align="v"), 
          nrow=1, common.legend=T, legend="none")

dev.off()

