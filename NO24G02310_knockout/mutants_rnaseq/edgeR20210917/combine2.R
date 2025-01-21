#!/usr/bin/env Rscript

all_genes <- read.table('../IMET1v2_nuclear_genes.txt',header=F,stringsAsFactors=F)$V1
M4_up <- read.table('M4_DEG_Up2.txt',header=T,stringsAsFactors=F)$DEG_Up
M4_down <- read.table('M4_DEG_Down2.txt',header=T,stringsAsFactors=F)$DEG_Down

#str(all_genes)
#str(M4_up)
#str(M4_down)

M4_DEG <- c()
for (g in all_genes) {
    deg <- "Not"
    if (g %in% M4_up) {
        deg <- "Up"
    } else if (g %in% M4_down) {
        deg <- "Down"
    }
    M4_DEG <- c(M4_DEG, deg)
}

M4_df <- data.frame(Gene=all_genes, DEG=M4_DEG)
#str(M4_df)

M6_up <- read.table('M6_DEG_Up2.txt',header=T,stringsAsFactors=F)$DEG_Up
M6_down <- read.table('M6_DEG_Down2.txt',header=T,stringsAsFactors=F)$DEG_Down

M6_DEG <- c()
for (g in all_genes) {
    deg <- "Not"
    if (g %in% M6_up) {
        deg <- "Up"
    } else if (g %in% M6_down) {
        deg <- "Down"
    }
    M6_DEG <- c(M6_DEG, deg)
}

M6_df <- data.frame(Gene=all_genes, DEG=M6_DEG)

write.table(M4_df, 'M4_DEGs2.txt', sep="\t", quote=F, row.names=F)
write.table(M6_df, 'M6_DEGs2.txt', sep="\t", quote=F, row.names=F)

