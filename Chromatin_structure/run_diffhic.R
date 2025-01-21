#!/usr/bin/env Rscript

library('diffHic')

frag <- cutGenome("IMET1v2.fasta","GATC",4)

param <- pairParam(frag)

diagnostics <- preparePairs("NannoH0_md_sort.bam", param, file="H0Rep1.h5", dedup=TRUE, minq=20)

names(diagnostics)
diagnostics$pairs
diagnostics$same.id
diagnostics$chimeras

diagnostics$chimeras[["invalid"]]/diagnostics$chimeras[["multi"]]

# find fragment length distribution
diags <- getPairData("H0Rep1.h5", param)
#hist(diags$length[diags$length < 1000], ylab="Frequency", xlab="Spacing (bp)", main="", col="grey80")

#min.inward <- 1000
#min.outward <- 25000
#prunePairs("H0Rep1.h5", param, file.out="H0Rep1_trimmed.h5", max.frag=600, min.inward=min.inward, min.outward=min.outward)

intra <- !is.na(diags$insert)
table(diags$orientation[!intra])

# Constructing the histograms for each orientation.
llinsert <- log2(diags$insert + 1L)
intra <- !is.na(llinsert)
breaks <- seq(min(llinsert[intra]), max(llinsert[intra]), length.out=30)
inward <- hist(llinsert[diags$orientation==1L], plot=FALSE, breaks=breaks)
outward <- hist(llinsert[diags$orientation==2L] ,plot=FALSE, breaks=breaks)
samestr <- hist(llinsert[diags$orientation==0L | diags$orientation==3L], plot=FALSE, breaks=breaks)
samestr$counts <- samestr$counts/2
# Setting up the axis limits.
ymax <- max(inward$counts, outward$counts, samestr$counts)/1e6
xmax <- max(inward$mids, outward$mids, samestr$mids)
xmin <- min(inward$mids, outward$mids, samestr$mids)
# Making a plot with all orientations overlaid.
plot(0,0,type="n", xlim=c(xmin, xmax), ylim=c(0, ymax), xlab=expression(log[2]~"[insert size (bp)]"), ylab="Frequency (millions)")
lines(inward$mids, inward$counts/1e6, col="darkgreen", lwd=2)
abline(v=log2(min.inward), col="darkgrey")
lines(outward$mids, outward$counts/1e6, col="red",lwd=2)
abline(v=log2(min.outward), col="darkgrey", lty=2)
lines(samestr$mids, samestr$counts/1e6, col="blue", lwd=2)
legend("topright", c("inward", "outward", "same"), col=c("darkgreen", "red", "blue"), lwd=2)

