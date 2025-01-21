#!/usr/bin/env Rscript

setwd('./')

set.seed(123456)
iter <- 1000

library(dplyr)
library(ggpubr)
library(extrafont)
#font_import()
loadfonts(device = "pdf")

data <- read.delim("H0_H24_CS_DEG.txt", header=T, stringsAsFactors=F)
d1 <- data[data$CS==1,]
d2 <- data[data$CS==2,]
d3 <- data[data$CS==3,]
d4 <- data[data$CS==4,]
d5 <- data[data$CS==5,]

dlist <- list("1"=d1, "2"=d2, "3"=d3, "4"=d4, "5"=d5)

cs_num <- data %>% count(Transition)

df <- data.frame(CS1=character(0),CS2=character(0),n=numeric(0),
                 Type=character(0),Percent=numeric(0),stringsAsFactors=F)

df_real <- data.frame(CS1=character(0),CS2=character(0),n=numeric(0),
                 Type=character(0),Percent=numeric(0),stringsAsFactors=F)

for (i in 1:nrow(cs_num)) {
  cst <- cs_num[i,"Transition"]
  csn <- cs_num[i,"n"]
  cs1 <- unlist(strsplit(cst,"_"))[1]
  cs2 <- unlist(strsplit(cst,"_"))[2]
  di <- dlist[[cs1]]
  cat("Transition: ", cst, "\n")
  if (csn >= 10) {
    ## actual percent of DE up and down
    datat <- data[data$Transition==cst,]
    up <- nrow(datat[datat$DEG == "Up",])*100.0/nrow(datat)
    dn <- nrow(datat[datat$DEG == "Down",])*100.0/nrow(datat)
    r1 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Up",Percent=up,stringsAsFactors=F)
    r2 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Down",Percent=dn,stringsAsFactors=F)
    df_real <- rbind(rbind(df_real, r1),r2)
    for (j in 1:iter) {
      sdj <- sample_n(di, csn)
      up_pct <- nrow(sdj[sdj$DEG == "Up",])*100.0/csn
      dn_pct <- nrow(sdj[sdj$DEG == "Down",])*100.0/csn
      item1 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Up",Percent=up_pct,stringsAsFactors=F)
      item2 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Down",Percent=dn_pct,stringsAsFactors=F)
      df <- rbind(rbind(df, item1),item2)
    }
  } else {
    r1 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Up",Percent=NA,stringsAsFactors=F)
    r2 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Down",Percent=NA,stringsAsFactors=F)
    df_real <- rbind(rbind(df_real, r1),r2)
    item1 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Up",Percent=NA,stringsAsFactors=F)
    item2 <- data.frame(CS1=paste0("CS",cs1),CS2=paste0("CS",cs2),n=csn,Type="Down",Percent=NA,stringsAsFactors=F)
    df <- rbind(rbind(df, item1),item2)
  }
}

##################################

pdf("NannoEpi_Fig5CS.pdf", useDingbats=FALSE, width=8, height=6, family="Arial")

df$Type <- factor(df$Type, levels=c('Up','Down'))
df$CS_LC <- df$CS2

ggboxplot(df, x="CS1", y="Percent", group="CS_LC", color="CS_LC", palette="jco",
              size=0.2, bxp.errorbar=T, facet.by="Type", ncol=1) + 
  geom_point(data=df_real, aes(shape=Type), fill="red", colour="red", size=2, 
             position=position_dodge2(width=0.75, preserve="single", padding=0)) +
  scale_shape_manual(values=c(24,25)) + theme_bw() + 
  theme(text=element_text(size=14, family = "Arial"), legend.position="right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("CS_HC") + ylab("Percent of DEGs (%)")

dev.off()

