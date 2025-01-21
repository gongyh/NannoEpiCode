library("ggvenn")

M4DEs <- read.table('M4_DEs.txt',header=F)$V1
M6DEs <- read.table('M6_DEs.txt',header=F)$V1

x <- list(M4=M4DEs, M6=M6DEs)

ggvenn(x,stroke_size=0.5)

