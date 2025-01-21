#!/usr/bin/env Rscript

library(HiCcompare)

# read in files
mat_hc1 <- read.table("../hic_results/matrix/H0rep/raw/rfbin/H0rep_rfbin.matrix")
bed_hc1 <- read.table("../hic_results/matrix/H0rep/raw/rfbin/H0rep_rfbin_abs.bed")
# convert to BEDPE
dat_hc1 <- hicpro2bedpe(mat_hc1, bed_hc1)
hc1 <- dat_hc1$cis # extract intrachromosomal matrices

# read in files
mat_hc2 <- read.table("../hic_results/matrix/H24rep/raw/rfbin/H24rep_rfbin.matrix")
bed_hc2 <- read.table("../hic_results/matrix/H24rep/raw/rfbin/H24rep_rfbin_abs.bed")
# convert to BEDPE
dat_hc2 <- hicpro2bedpe(mat_hc2, bed_hc2)
hc2 <- dat_hc2$cis # extract intrachromosomal matrices

# for all chromosomes
hic.list <- mapply(create.hic.table, hc1, hc2, SIMPLIFY = FALSE, scale = FALSE)

# Total Sum Scaling
hic.list <- total_sum(hic.list)

# Joint Normalization
# Multiple hic.tables can be processed in parallel by entering a list of hic.tables
hic.list <- hic_loess(hic.list, parallel = TRUE)

# Difference Detection
hic.table <- list()
for (i in 1:30) {
  hic.table[[i]] <- hic_compare(hic.list[[i]], adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)
  num <- dim(hic.table[[i]][hic.table[[i]]$p.adj<0.01,])[1]
  str(num)
}

hict <- as.data.frame(do.call(rbind, hic.table))
str(hict[hict$p.adj<0.01,])

write.table(hict[hict$p.adj<0.01,], file="Rep2_rfbin.txt", quote=F, sep="\t")

