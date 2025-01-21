#!/usr/bin/env Rscript

library(HiCcompare)
library(BiocParallel)

# read in files
mat_hc1 <- read.table("../hic_results/matrix/H0rep/raw/1000/H0rep_1000.matrix")
bed_hc1 <- read.table("../hic_results/matrix/H0rep/raw/1000/H0rep_1000_abs.bed")
# convert to BEDPE
dat_hc1 <- hicpro2bedpe(mat_hc1, bed_hc1)
hc1 <- dat_hc1$cis # extract intrachromosomal matrices

# read in files
mat_hc2 <- read.table("../hic_results/matrix/H24rep/raw/1000/H24rep_1000.matrix")
bed_hc2 <- read.table("../hic_results/matrix/H24rep/raw/1000/H24rep_1000_abs.bed")
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
hic.table <- hic_compare(hic.list, adjust.dist = TRUE, p.method = 'fdr', parallel=T, BP_param = MulticoreParam(workers=30))

hict <- as.data.frame(do.call(rbind, hic.table))
str(hict[hict$p.adj<0.01,])

write.table(hict[hict$p.adj<0.01,], file="Rep2_1k.txt", quote=F, sep="\t")

