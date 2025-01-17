#!/usr/bin/env Rscript

fs <- list.files(getwd(), pattern=".lc.txt")

nsamples <- length(fs)

lc <- list()
max_x <- 0
max_y <- 0
if (length(fs)>0) {
    i <- 1
    for (f in fs) {
        f_lc <- read.delim(f)
        f_lc_mx <- max(f_lc$TOTAL_READS)
	f_lc_my <- max(f_lc$EXPECTED_DISTINCT)
	if (f_lc_mx>max_x) max_x<-f_lc_mx
	if (f_lc_my>max_y) max_y<-f_lc_my
        lc[[i]] <- f_lc
	i <- i+1
    }
}

pdf("preseq.pdf")

colors <- rainbow(nsamples)

plot(range(0,max_x),range(0,max_y),type='n',xlab='TOTAL_READS',ylab='EXPECTED_DISTINCT')

for (i in 1:nsamples) {
    lines(lc[[i]]$TOTAL_READS,lc[[i]]$EXPECTED_DISTINCT,col=colors[i],type='l',lwd=1.5)
}

legend(0.7*max_x, 0.3*max_y, legend=fs, cex=0.8, col=colors, lwd=1.5)

dev.off()

