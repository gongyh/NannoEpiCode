#!/usr/bin/env Rscript

setwd('./')

set.seed(123456)
iter <- 1000

library(dplyr)
library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

library(gtools)

data <- read.delim("H0_H24_CS_DEG.txt", header=T, stringsAsFactors=F)
d1 <- data[data$CS==1,]
d2 <- data[data$CS==2,]
d3 <- data[data$CS==3,]
d4 <- data[data$CS==4,]
d5 <- data[data$CS==5,]

dlist <- list("1"=d1, "2"=d2, "3"=d3, "4"=d4, "5"=d5)

cs_num <- data %>% count(Transition)

df_real <- data.frame(CS1=character(0),CS2=character(0),n=numeric(0),Type=character(0),
  Percent=numeric(0),HC_percent=numeric(0),pvalue=numeric(0),stringsAsFactors=F)

for (i in 1:nrow(cs_num)) {
  cst <- cs_num[i,"Transition"]
  csn <- cs_num[i,"n"]
  cs1 <- unlist(strsplit(cst,"_"))[1]
  cs2 <- unlist(strsplit(cst,"_"))[2]
  di <- dlist[[cs1]]
  cat("Transition: ", cst, "\n")
  if ( (cs1!=cs2) && (csn >= nrow(di)*0.05) ) {
    ## actual percent of DE up and down
    datat <- data[data$Transition==cst,]
    up <- nrow(datat[datat$DEG == "Up",])*100.0/nrow(datat)
    dn <- nrow(datat[datat$DEG == "Down",])*100.0/nrow(datat)
    data_HC <- data[data$CS==cs1,]
    HC_up <- nrow(data_HC[data_HC$DEG == "Up",])*100.0/nrow(data_HC)
    HC_dn <- nrow(data_HC[data_HC$DEG == "Down",])*100.0/nrow(data_HC)
    up_all <- c()
    dn_all <- c()
    for (j in 1:iter) {
      sdj <- sample_n(di, csn)
      up_pct <- nrow(sdj[sdj$DEG == "Up",])*100.0/csn
      dn_pct <- nrow(sdj[sdj$DEG == "Down",])*100.0/csn
      up_all <- c(up_all, up_pct)
      dn_all <- c(dn_all, dn_pct)
    }
    up_pvalue <- min(sum(up_all>up)/iter, sum(up_all<up)/iter)
    dn_pvalue <- min(sum(dn_all>dn)/iter, sum(dn_all<dn)/iter)
    r1 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Up",
            Percent=up,HC_percent=HC_up,pvalue=up_pvalue,stringsAsFactors=F)
    r2 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Down",
            Percent=dn,HC_percent=HC_dn,pvalue=dn_pvalue,stringsAsFactors=F)
    df_real <- rbind(rbind(df_real, r1),r2)
  }
}


##################################

pdf("NannoEpi_Fig5CS_bar.pdf", useDingbats=FALSE, width=4, height=7, family="Arial")

df_real$CS1 <- factor(df_real$CS1, levels=c('CS1','CS2','CS3','CS4','CS5'), ordered=T)
df_real$CS2 <- factor(df_real$CS2, levels=c('CS1','CS2','CS3','CS4','CS5'), ordered=T)
df_real$Type <- factor(df_real$Type, levels=c('Down','Up'), ordered=T)
df_real$CS_LC <- df_real$CS2
df_real <- arrange(df_real, CS1, CS_LC)
df_real$CS_CS <- paste0(df_real$CS1,"_",df_real$CS2)
df_real$CS_CS <- factor(df_real$CS_CS, ordered=T,
    levels=c('CS2_CS1','CS3_CS1','CS1_CS2','CS4_CS3','CS3_CS4','CS5_CS4','CS4_CS5'))

str(df_real)

ggbarplot(df_real, x="Type", y="Percent", fill="Type", color="Type", palette="jco", 
          facet.by="CS_CS", ncol=1, strip.position="left", width=0.5) + 
  scale_y_continuous(position="right", limits=c(0,35)) +
  geom_point(aes(y=HC_percent), fill="red", colour="red", size=2) +
  geom_text(aes(label = stars.pval(pvalue)), vjust=0.75, hjust=-1) +
  theme_bw() + coord_flip() +
  theme(text=element_text(size=14, family = "Arial"), legend.position="bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(color="black")) +
  xlab("") + ylab("Percent of DEGs (%)")

dev.off()

