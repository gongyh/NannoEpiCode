#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

HC_R1_df <- read.table('HC_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)
LC_R1_df <- read.table('LC_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)

bin_gene <- read.table('bins_40000_genes.bed', sep="\t", header=F)
gene_log2fc <- read.table('LC_logfc.tsv', sep="\t", header=F)
gene_deg <- read.table('H0_H24_DEG.txt', sep="\t", header=T)
gene_tpm <- read.table('HC_VLC.TMM.TPM.matrix.noCpMt',sep="\t", header=T)

bin_gene_tpm <- merge(bin_gene,gene_log2fc,by.x="V4",by.y="V1")
#str(bin_gene_log2fc)
bin_gene_log2fc <- merge(bin_gene_tpm, gene_tpm, by.x="V4", by.y="Genes")
bin_gene_log2fc$pos <- (bin_gene_log2fc$V2.x + bin_gene_log2fc$V3)/2

bin_gene_log2fc <- merge(bin_gene_log2fc,gene_deg,by.x="V4",by.y="Gene")
str(bin_gene_log2fc)

HC_R1_df$pos <- (HC_R1_df$V2 + HC_R1_df$V3)/2
LC_R1_df$pos <- (LC_R1_df$V2 + LC_R1_df$V3)/2

HC_R1_df$sign <- HC_R1_df$V4 >= 0.0
LC_R1_df$sign <- LC_R1_df$V4 >= 0.0

pdf("NannoEpi_HiC_AB_HCLC_log2fc.pdf", useDingbats=FALSE, width=6, height=5, family="Arial", onefile=F)

HC1_data <- HC_R1_df
LC1_data <- LC_R1_df

bgl_data <- bin_gene_log2fc
#str(bgl_data)

HC1_data$grp <- 'HC'
LC1_data$grp <- 'LC'

data <- rbind(HC1_data, LC1_data)
data$chr <- as.numeric(sub("chr","",data$V1))
data$sign <- factor(data$sign, levels=c(TRUE,FALSE),ordered=T)
str(data)

HCLC_data <- HC1_data
HCLC_data$grp <- 'HC*LC'
HCLC_data$V4 <- HC1_data$V4*LC1_data$V4
HCLC_data$sign <- HCLC_data$V4 >= 0.0
HCLC_data$chr <- as.numeric(sub("chr","",HCLC_data$V1))
HCLC_data$sign <- factor(HCLC_data$sign, levels=c(TRUE,FALSE),ordered=T)
dataDiff <- rbind(data,HCLC_data)
dataDiff$grp <- factor(dataDiff$grp, levels=c("HC","LC","HC*LC"), ordered=T)

dataDiff <- dataDiff[dataDiff$grp %in% c("HC","LC"),]

a <- ggdotchart(dataDiff, x="pos", y="V4", sorting="none", add = "segment", add.params=list(size=0.5,color="sign"),
           size=0.05, color="sign", dot.size=0.01, palette=c("red","blue")) + 
     facet_grid(vars(grp),vars(chr),scales="free_x",space="free_x",switch="x") +
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), panel.spacing.x = unit(0, "cm"), 
        strip.text.x = element_text(angle=90,size=6), strip.background=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.y = element_text(color="black")) +
  xlab("") + ylab("E1")

print("a")

bgl_data$HC <- bgl_data$Control
bgl_data$LC <- bgl_data$VLC
ab24 <- melt(bgl_data,measure.vars=c("HC","LC"),id.vars=c("V1","pos"),variable.name="grp")
str(ab24)
ab24_df <- merge(ab24,data,by=c("V1","pos","grp"))
str(ab24_df)

cmps<-list(c("TRUE","FALSE"))

ab <- ggbarplot(ab24_df, x="sign", y="value", group="sign", fill="sign", palette=c("red","blue"), facet.by="grp",nrow=1,
                add="mean_se", error.plot="errorbar", outlier.shape=NA) + coord_cartesian(ylim = c(0, 250)) +
  scale_x_discrete(labels=c("TRUE"="A","FALSE"="B")) + theme_bw() +
  theme(text=element_text(size=9, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Compartment") + ylab("RNA-Seq (TPM)") + 
  stat_compare_means(comparisons=cmps, method="wilcox.test", aes(labels=..p.signif..),
                     label.y=200, tip.length=0.003)

print("ab")

data2 <- merge(HC1_data,LC1_data,by=c("V1","pos"))
df <- data2[,c("V1","pos","sign.x","sign.y","V4.x","V4.y")]
df$type <- paste(data2$sign.x,data2$sign.y,sep="_")
str(df)
df <- df[df$V4.x*df$V4.y != 0,]
str(df)
str(bgl_data)
data3 <- merge(bgl_data,df)
data3$type <- factor(data3$type, levels=c("TRUE_FALSE","TRUE_TRUE","FALSE_TRUE","FALSE_FALSE"), ordered=T)
cmps2 <- list(c("TRUE_FALSE","TRUE_TRUE"),c("FALSE_TRUE","FALSE_FALSE"))

chosed <- c("chr1","chr2","chr3","chr4","chr5","chr9","chr10",
            "chr11","chr12","chr15","chr16","chr17","chr18","chr19",
            "chr23","chr24","chr25","chr27","chr28","chr30")
data35 <- data3[data3$V1 %in% chosed,]

str(data35)

b <- ggviolin(data35, x="type", y="V2.y", fill="type", palette="Paired", size=0.4, width=0.7) + 
     geom_boxplot(color="black", fill = "white", width=0.2, outlier.size=0.3) +
  scale_x_discrete(labels=c("TRUE_FALSE"="A2B","FALSE_TRUE"="B2A", "TRUE_TRUE"="A2A", "FALSE_FALSE"="B2B")) +
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Compartment changes") + ylab("RNA-Seq (log2fc)") +
  stat_compare_means(comparisons=cmps2, method="wilcox.test",aes(labels=..p.signif..))


data35$chr <- as.numeric(sub("chr","",data35$V1))
e <- ggboxplot(data35, x="type", y="V2.y", notch=T, outlier.shape=NA, bxp.errorbar=T,facet.by="chr", ncol=8) +
  scale_x_discrete(labels=c("TRUE_FALSE"="A2B","FALSE_TRUE"="B2A", "TRUE_TRUE"="A2A", "FALSE_FALSE"="B2B")) +
  theme_bw() +
  theme(text=element_text(size=8, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("RNA-Seq (log2fc)") +
  stat_compare_means(comparisons=cmps2, method="t.test",size=2)

ggarrange(a,ggarrange(ggplot(),ab,nrow=1,widths=c(1.75,1)), ncol=1)

dev.off()


pdf("NannoEpi_HiC_A2B_HCLC_TPM.pdf", useDingbats=FALSE, width=2.25, height=5, family="Arial", onefile=F)
d35 <- data35
dpct <- d35 %>% count(type, DEG)

d35$cs_diff <- d35$V4.y-d35$V4.x
cmps3 <- list(c("Up","Not"), c("Down","Not"))
d35 %>% group_by(DEG) %>% summarise_at(vars(cs_diff), list(name = mean))
p1 <- ggviolin(d35, x="DEG", y="cs_diff", add="boxplot", outlier.shape=NA, bxp.errorbar=T) +
  theme_bw() + theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("RNA-Seq changes") + ylab("Compartment score changes") +
  stat_compare_means(comparisons=cmps3, method="wilcox.test",aes(labels=..p.signif..))

abud <- data.frame(type=c("A2B","A2A","B2A","B2B"),ratio=c(7/10,636/475,65/35,584/423))
abud$type <- factor(abud$type, levels=c("A2B","A2A","B2A","B2B"), ordered=T)
pud <- ggbarplot(abud,x="type",y="ratio",fill="type", palette="Paired") +
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("Compartment changes") + ylab("#DEG_Up / #DEG_Down")

ggarrange(pud,b,ncol=1)

dev.off()


