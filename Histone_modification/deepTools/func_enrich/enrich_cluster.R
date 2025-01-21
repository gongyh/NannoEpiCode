#!/usr/bin/env Rscript

library("optparse")
 
option_list = list(
  make_option(c("-g", "--group"), type="character", default=NULL, 
              help="H0 or H24", metavar="character"),
  make_option(c("-m", "--mark"), type="character", default=NULL,
              help="Histone mark: H3k9ac, H3k27ac, H3kcr, H3k4m2", metavar="character"),
  make_option(c("-c", "--cluster"), type="character", default=NULL,
              help="cluster_1, cluster_2, cluster_3, cluster_4, cluster_5", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#H24_H3kcr_k5.plotHeatmap.sorted.bed
bed_file <- paste0("../",opt$group,"_",opt$mark,"_k5.plotHeatmap.sorted.bed")

if (file.exists(bed_file)){
  print("BEGIN ......")
} else {
  print_help(opt_parser)
  stop("No heatmap bed file found. Wrong parameter?", call.=FALSE)
}

data <- read.delim(bed_file, stringsAsFactors=F)
#str(data)

library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(DOSE)
library(dplyr)
library(org.Noceanica.eg.db)

genes <- data[data$deepTools_group==opt$cluster,"name"]
genes <- stringr::str_replace_all(genes,".1","")
str(genes)

res_bp <- enrichGO(genes, 'org.Noceanica.eg.db', ont="BP", pvalueCutoff=0.01, keyType = "GID")
res_bp2 <- simplify(res_bp)
res_mf <- enrichGO(genes, 'org.Noceanica.eg.db', ont="MF", pvalueCutoff=0.01, keyType = "GID")
res_mf2 <- simplify(res_mf)
res_cc <- enrichGO(genes, 'org.Noceanica.eg.db', ont="CC", pvalueCutoff=0.01, keyType = "GID")
res_cc2 <- simplify(res_cc)

res_table <- rbind(res_bp2@result, res_mf2@result, res_cc2@result)
data <- res_table[res_table$p.adjust<0.01,]
data$Ontology <- Ontology(GOTERM[data$ID])

str(data)

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

pdf(opt$out, useDingbats = FALSE, width = 10, height = 5, family="Arial")

p1 <- ggbarplot(data, x="ID", y="qvalue", fill="Ontology", color="Ontology", group="Ontology", sort.val="desc",
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

df <- data
id <- df$ID
des <- df$Description
glist <- apply(data,1,function(x){strsplit(x["geneID"],"/")[[1]]},simplify=F)
#names(glist) <- des
d <- list2df(glist)
res <- tibble::tibble(Description = split(d[,1], d[,2]))
p2 <- ggplot(res, aes_(x = ~Description)) + geom_bar() +
    geom_text(stat='count', aes(label=after_stat(count)), vjust=1, colour="white") +
    theme_dose(font.size = 12) +
    xlab(NULL) + ylab(NULL) +
    ggupset::scale_x_upset(order_by = "freq")

ggarrange(p1,p2, nrow=1)

dev.off()

