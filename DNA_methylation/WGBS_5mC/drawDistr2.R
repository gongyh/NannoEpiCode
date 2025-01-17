#!/usr/bin/env Rscript

setwd('./')

#setwd(dirname(sys.frame(1)$ofile))
#setwd('E:\\工作相关\\NannoGenome\\5mC')

df <- t(read.delim('Hdistr_fix2.txt',row.names=1, header=FALSE))

ss_df <- df[,c("ss_None","ss_Low","ss_Medium","ss_High")]
#nss_df <- df[,c("nss_None","nss_Low","nss_Medium","nss_High")]

df2 <- t(read.delim('Ldistr_fix2.txt',row.names=1, header=FALSE))

ss_df2 <- df2[,c("ss_None","ss_Low","ss_Medium","ss_High")]
#nss_df2 <- df2[,c("nss_None","nss_Low","nss_Medium","nss_High")]

#par(family="Arial",ps=20)

library(extrafont)
loadfonts(device = "pdf")

pdf("MLdistr_fix2.pdf", width=10, height=4, family="Arial", useDingbats = FALSE)

par(mfcol=c(1,2), mar=c(2,4,1,1))

plot(1:50,ss_df[,4],col="red",type='l',ylim=c(0.05,0.25),xaxt="n",xlab="",ylab="5mC level (%)")
points(1:50,ss_df[,3],col='green',type='l')
points(1:50,ss_df[,2],col='blue',type='l')
points(1:50,ss_df[,1],col='gray',type='l')
abline(v=c(15,36),col="black",lty='dashed')
axis(1,at=c(1,15,36,50),labels=c("1.5kb","TSS","TES","1.5kb"))
legend(17,0.25,legend=c('High','Medium','Low','Silent'),col=c("red",'green','blue','gray'),lty=1,
       bty="n", title="HC")

#plot(1:50,nss_df[,4],col="red",type='l',ylim=c(0.05,0.25),xaxt="n",xlab="",ylab="5mC level (%)")
#points(1:50,nss_df[,3],col='green',type='l')
#points(1:50,nss_df[,2],col='blue',type='l')
#points(1:50,nss_df[,1],col='gray',type='l')
#abline(v=c(15,36),col="black",lty='dashed')
#axis(1,at=c(1,15,36,50),labels=c("1.5kb","TSS","TES","1.5kb"))
#legend(1,0.25,legend=c('High','Medium','Low','None'),col=c("red",'green','blue','gray'),lty=1)

plot(1:50,ss_df2[,4],col="red",type='l',ylim=c(0.05,0.25),xaxt="n",xlab="",ylab="5mC level (%)")
points(1:50,ss_df2[,3],col='green',type='l')
points(1:50,ss_df2[,2],col='blue',type='l')
points(1:50,ss_df2[,1],col='gray',type='l')
abline(v=c(15,36),col="black",lty='dashed')
axis(1,at=c(1,15,36,50),labels=c("1.5kb","TSS","TES","1.5kb"))
legend(17,0.25,legend=c('High','Medium','Low','Silent'),col=c("red",'green','blue','gray'),lty=1,
       bty="n", title="LC")

#plot(1:50,nss_df2[,4],col="red",type='l',ylim=c(0.05,0.25),xaxt="n",xlab="",ylab="5mC level (%)")
#points(1:50,nss_df2[,3],col='green',type='l')
#points(1:50,nss_df2[,2],col='blue',type='l')
#points(1:50,nss_df2[,1],col='gray',type='l')
#abline(v=c(15,36),col="black",lty='dashed')
#axis(1,at=c(1,15,36,50),labels=c("1.5kb","TSS","TES","1.5kb"))
#legend(1,0.25,legend=c('High','Medium','Low','None'),col=c("red",'green','blue','gray'),lty=1)

dev.off()
