#!/bin/env Rscript

source("TopDom_v0.0.2.R")

#~/bin/HiC-Pro_2.9.0/bin/utils/sparseToDense.py -c -d -b matrix/NannoH0/raw/5000/NannoH0_5000_abs.bed -o test.matrix matrix/NannoH0/iced/5000/NannoH0_5000_iced.matrix

TopDom(matrix.file="chr1_test.matrix",window.size=10,outFile="test")

