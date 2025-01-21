#!/bin/env Rscript

#################################################################
# Function: Compare compartment and TAD based on HiC contact matrix
# Call: Rscript compare.R -m1 1.matrix -m2 2.matrix -b bed_file -s chr1_NannoH0_5000.is30001.ids20001.insulation -S chr1_NannoH24_5000.is30001.ids20001.insulation [-c chr1 -l1 Sample1 -l2 Sample2 -o outfile]
# R packages used: optparse, HiTC, lattice, reshape2
# Email: Yanhai Gong, gongyh@qibebt.ac.cn
#################################################################

## install necessary libraries
p <- c("optparse","lattice","reshape2")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

#install.packages("devtools")
#devtools::install_bitbucket("graumannlabtools/multipanelfigure")

## clean R environment
rm(list = ls())
setwd('./')

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)

# make option list and parse command line
option_list <- list(
  make_option(c("-i", "--matrix_file1"), type="character", help="matrix 1 [Required]"),
  make_option(c("-I", "--matrix_file2"), type="character", help="matrix 2 [Required]"),
  make_option(c("-s", "--insulation_file1"), type="character", help="insulation file1 [Required]"),
  make_option(c("-S", "--insulation_file2"), type="character", help="insulation file2 [Required]"),
  make_option(c("-b", "--bed_file"), type="character", help="bed file [Required]"),
  make_option(c("-r", "--resolution"), type="integer", help="resolution [Required]"),
  make_option(c("-c", "--chr_name"), type="character", default='chr1', help="chromosome name [Optional, default %default]"),
  make_option(c("-o", "--outfile"), type="character", default='out.pdf', help="Output figure [Optional, default %default]"),
  make_option(c("-l", "--label1"), type="character", default='Sample1', help="label1 [Optional, default %default]"),
  make_option(c("-L", "--label2"), type="character", default='Sample2', help="label2 [Optional, default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$matrix_file1)) stop('Please input a matrix file')
if(is.null(opts$matrix_file2)) stop('Please input another matrix file')
if(is.null(opts$bed_file)) stop('Please input a bed file')

m1f <- opts$matrix_file1
m2f <- opts$matrix_file2
is1f <- opts$insulation_file1
is2f <- opts$insulation_file2
bedf <- opts$bed_file
chr_name <- opts$chr_name
outf <- opts$outfile
label1 <- opts$label1
label2 <- opts$label2

bed_df <- read.table(bedf,header=F,sep="\t")
bed_b_e <- bed_df[bed_df$V1==chr_name,2:3]
starts <- bed_b_e$V2
chr_length <- dplyr::last(bed_b_e$V3)
resolution <- opts$resolution
x<-rep(starts,each=resolution)[0:chr_length]

library('HiTC')

hiC1 <- importC(con=m1f, xgi.bed=bedf, ygi.bed=bedf,rm.trans=TRUE)
hiC2 <- importC(con=m2f, xgi.bed=bedf, ygi.bed=bedf,rm.trans=TRUE)

hiC1_chr <- hiC1[[paste(chr_name,chr_name,sep="")]]
hiC2_chr <- hiC2[[paste(chr_name,chr_name,sep="")]]

class(hiC1_chr@intdata)<-"dsCMatrix"
hiC1_chr@intdata@uplo <- "U"

class(hiC2_chr@intdata)<-"dsCMatrix"
hiC2_chr@intdata@uplo <- "U"

#library(multipanelfigure)

pdf(outf)
opar <- par() # store par settings

#figure1 <- multi_panel_figure(
#  width = 18, height = 36,
#  columns = 2, rows = 6)

#a_base_plot <- capture_base_plot( mapC(hiC1_chr, title=paste(label1,chr_name,sep=" - ")) )
#b_base_plot <- capture_base_plot( mapC(hiC2_chr, title=paste(label2,chr_name,sep=" - ")) )

#figure1 %<>% fill_panel(a_base_plot, column = 1)
#figure1 %<>% fill_panel(a_base_plot, column = 2)
#figure1 %<>% fill_panel(b_base_plot, column = 1)
#figure1 %<>% fill_panel(b_base_plot, column = 2)

mapC(hiC1_chr, hiC2_chr, title=paste(paste(label1,label2,sep=" vs "),chr_name,sep=" , "))

#par(mfrow=c(2,1))

pc1 <- pca.hic(hiC1_chr, normPerExpected=TRUE, npc=1)
#plot(start(pc1$PC1), score(pc1$PC1), type="h", xlab=paste(label1,chr_name,sep=" , "), ylab="PC1vec", frame=FALSE)

pc2 <- pca.hic(hiC2_chr, normPerExpected=TRUE, npc=1)
#plot(start(pc2$PC1), score(pc2$PC1), type="h", xlab=paste(label2,chr_name,sep=" , "), ylab="PC1vec", frame=FALSE)

#par(opar)

library(lattice)
library(reshape2)

compartment1 <- data.frame(loci=start(pc1$PC1),Sample1=score(pc1$PC1))
compartment2 <- data.frame(loci=start(pc2$PC1),Sample2=score(pc2$PC1))
compartment <- merge(compartment1,compartment2,by="loci")
mm <- melt(compartment,id.var="loci",value.name="PC1")

levels(mm$variable)[levels(mm$variable)=="Sample1"] <- label1
levels(mm$variable)[levels(mm$variable)=="Sample2"] <- label2

mm$variable <- factor(mm$variable,levels=c(label2,label1),ordered=TRUE)

cp <- xyplot(PC1~loci|variable,data=mm,type="h",scales=list(y=list(relation="free")),layout=c(1,2),
  xlab=chr_name,grid=TRUE)

cp

#figure1 %<>% fill_panel(cp,column=1,row=c(3,4))

###########
## draw RNA-seq coverage hist
par(opar)
h0_1_df <- read.table("../../0h-1.cov",header=F,sep="\t")
h0_1_chr_df <- h0_1_df[h0_1_df$V1==chr_name,c(1,3)]
h0_1_chr_df$x <- x
h0_1_chr_hist <- tapply(h0_1_chr_df$V3,h0_1_chr_df$x,sum)/resolution

h0_2_df <- read.table("../../0h-2.cov",header=F,sep="\t")
h0_2_chr_df <- h0_2_df[h0_2_df$V1==chr_name,c(1,3)]
h0_2_chr_df$x <- x
h0_2_chr_hist <- tapply(h0_2_chr_df$V3,h0_2_chr_df$x,sum)/resolution

h24_1_df <- read.table("../../VLC24-1.cov",header=F,sep="\t")
h24_1_chr_df <- h24_1_df[h24_1_df$V1==chr_name,c(1,3)]
h24_1_chr_df$x <- x
h24_1_chr_hist <- tapply(h24_1_chr_df$V3,h24_1_chr_df$x,sum)/resolution

h24_2_df <- read.table("../../VLC24-2.cov",header=F,sep="\t")
h24_2_chr_df <- h24_2_df[h24_2_df$V1==chr_name,c(1,3)]
h24_2_chr_df$x <- x
h24_2_chr_hist <- tapply(h24_2_chr_df$V3,h24_2_chr_df$x,sum)/resolution

par(mfrow=c(4,1))
par(mar=c(5,4,0.5,2))
par(oma=c(0,0,0,0))
plot(starts,h0_1_chr_hist,type="h",xlab="loci (bp)",ylab="H0_rep1")
plot(starts,h0_2_chr_hist,type="h",xlab="loci (bp)",ylab="H0_rep2")
plot(starts,h24_1_chr_hist,type="h",xlab="loci (bp)",ylab="H24_rep1")
plot(starts,h24_2_chr_hist,type="h",xlab="loci (bp)",ylab="H24_rep2")
###########

par(opar) # restore default, to avoid influence from lattice xyplot

# read insulation files
is1_raw <- read.delim(is1f)
is2_raw <- read.delim(is2f)

loci <-is1_raw$midpoint
is1 <- is1_raw$insulationScore
is2 <- is2_raw$insulationScore
isDiff <- is1_raw$insulationScore-is2_raw$insulationScore

par(mfrow=c(2,1))

#is1_base_plot <- capture_base_plot({
  plot(loci[!is.na(is1)],is1[!is.na(is1)],type="l",col="red",ylab="insulation score",xlab=chr_name,ylim=c(-1.5,max(is1[!is.na(is1)])+0.1))
  lines(loci[!is.na(is2)],is2[!is.na(is2)],type="l",col="blue")
  abline(h=0)
#})

#is2_base_plot <- capture_base_plot({
  plot(loci[!is.na(isDiff)],isDiff[!is.na(isDiff)],type="h",col="grey",xlab=chr_name,ylab="difference of IS",ylim=c(-0.5,0.5))
  abline(h=0)
#})

#di1<-directionalityIndex(hiC1_chr,winup=25e3,windown=25e3)
#di2<-directionalityIndex(hiC2_chr,winup=25e3,windown=25e3)

#par(mfrow=c(2,1))

#di1_base_plot <- capture_base_plot({
  #barplot(di1, col=ifelse(di1>0,"darkred","darkgreen"),main=paste("DI of ",label1,sep=""),ylim=c(-100,100)) 
#})

#di2_base_plot <- capture_base_plot({
  #barplot(di2, col=ifelse(di2>0,"darkred","darkgreen"),main=paste("DI of ",label2,sep=""),ylim=c(-100,100)) 
#})

#figure1 %<>% fill_panel(is1_base_plot, column = 2)
#figure1 %<>% fill_panel(is2_base_plot, column = 2)
#figure1 %<>% fill_panel(di1_base_plot, column = 2)
#figure1 %<>% fill_panel(di2_base_plot, column = 2)

#figure1

dev.off()

q()

