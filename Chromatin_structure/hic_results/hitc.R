#!/bin/env Rscript

library('HiTC')

hiC <- importC(con="matrix/NannoH0/iced/5000/NannoH0_5000_iced.matrix",xgi.bed="matrix/NannoH0/raw/5000/NannoH0_5000_abs.bed",
  ygi.bed="matrix/NannoH0/raw/5000/NannoH0_5000_abs.bed",rm.trans=FALSE)


#pdf('vis.pdf')
##mapC(hiC)
#mapC(hiC$chr1chr1)
#mapC(hiC$chr2chr2)
#mapC(hiC$chr3chr3)
#dev.off()

#pdf('cqc5k.pdf')
#par(mfrow=c(2,2))
#CQC(hiC, winsize = 1e+04, dev.new=FALSE, hist.dist=FALSE)
#dev.off()

pdf('pca5k.pdf')
pc <- pca.hic(hiC$chr1chr1, normPerExpected=FALSE, npc=1)
plot(start(pc$PC1), score(pc$PC1), type="h", xlab="chr1", ylab="PC1vec", frame=FALSE)
pc <- pca.hic(hiC$chr2chr2, normPerExpected=FALSE, npc=1)
plot(start(pc$PC1), score(pc$PC1), type="h", xlab="chr2", ylab="PC1vec", frame=FALSE)
dev.off()

