#!/usr/bin/env Rscript

lc_0_1 <- read.delim("IMET1v2.0_1.lc.txt")
lc_0_2 <- read.delim("IMET1v2.0_2.lc.txt")
lc_0_3 <- read.delim("IMET1v2.0_3.lc.txt")

lc_24_1 <- read.delim("IMET1v2.24_1.lc.txt")
lc_24_2 <- read.delim("IMET1v2.24_2.lc.txt")
lc_24_3 <- read.delim("IMET1v2.24_3.lc.txt")

max_x <- max(c(lc_0_1$TOTAL_READS[length(lc_0_1$TOTAL_READS)],lc_0_2$TOTAL_READS[length(lc_0_2$TOTAL_READS)],
	   lc_0_3$TOTAL_READS[length(lc_0_3$TOTAL_READS)],lc_24_1$TOTAL_READS[length(lc_24_1$TOTAL_READS)],
	   lc_24_2$TOTAL_READS[length(lc_24_2$TOTAL_READS)],lc_24_3$TOTAL_READS[length(lc_24_3$TOTAL_READS)]))

max_y <- max(c(lc_0_1$EXPECTED_DISTINCT,lc_0_2$EXPECTED_DISTINCT,lc_0_3$EXPECTED_DISTINCT,
	       lc_24_1$EXPECTED_DISTINCT,lc_24_2$EXPECTED_DISTINCT,lc_24_3$EXPECTED_DISTINCT))

pdf("lc.pdf")

colors <- rainbow(6)

plot(range(0,max_x),range(0,max_y),type='n',xlab='TOTAL_READS',ylab='EXPECTED_DISTINCT')
lines(lc_0_1$TOTAL_READS,lc_0_1$EXPECTED_DISTINCT,col=colors[1],type='l',lwd=1.5)
lines(lc_0_2$TOTAL_READS,lc_0_2$EXPECTED_DISTINCT,col=colors[2],type='l',lwd=1.5)
lines(lc_0_3$TOTAL_READS,lc_0_3$EXPECTED_DISTINCT,col=colors[3],type='l',lwd=1.5)
lines(lc_24_1$TOTAL_READS,lc_24_1$EXPECTED_DISTINCT,col=colors[4],type='l',lwd=1.5)
lines(lc_24_2$TOTAL_READS,lc_24_2$EXPECTED_DISTINCT,col=colors[5],type='l',lwd=1.5)
lines(lc_24_3$TOTAL_READS,lc_24_3$EXPECTED_DISTINCT,col=colors[6],type='l',lwd=1.5)
legend(0.8*max_x, 0.3*max_y, legend=c("lc_0_1","lc_0_2","lc_0_3","lc_24_1","lc_24_2","lc_24_3"),
       cex=0.8, col=colors, lwd=1.5)
dev.off()

