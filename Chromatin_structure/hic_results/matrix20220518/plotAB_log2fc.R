#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)
library(ggpubr)

library(extrafont)
#font_import()
loadfonts(device = "pdf")

HC_R1_df <- read.table('NannoH0_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)
HC_R2_df <- read.table('H0rep_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)
LC_R1_df <- read.table('NannoH24_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)
LC_R2_df <- read.table('H24rep_40000.cis.E1.bedgraph', sep="\t", header=F) %>% replace(is.na(.), 0)

bin_gene <- read.table('bins_40000_genes.bed', sep="\t", header=F)
gene_log2fc <- read.table('LC_logfc.tsv', sep="\t", header=F)

gene_tpm <- read.table('HC_VLC.TMM.TPM.matrix.noCpMt',sep="\t", header=T)

bin_gene_tpm <- merge(bin_gene,gene_log2fc,by.x="V4",by.y="V1")
#str(bin_gene_log2fc)
bin_gene_log2fc <- merge(bin_gene_tpm, gene_tpm, by.x="V4", by.y="Genes")
bin_gene_log2fc$pos <- (bin_gene_log2fc$V2.x + bin_gene_log2fc$V3)/2

HC_R1_df$pos <- (HC_R1_df$V2 + HC_R1_df$V3)/2
HC_R2_df$pos <- (HC_R2_df$V2 + HC_R2_df$V3)/2
LC_R1_df$pos <- (LC_R1_df$V2 + LC_R1_df$V3)/2
LC_R2_df$pos <- (LC_R2_df$V2 + LC_R2_df$V3)/2

HC_R1_df$sign <- HC_R1_df$V4 >= 0.0
HC_R2_df$sign <- HC_R2_df$V4 >= 0.0
LC_R1_df$sign <- LC_R1_df$V4 >= 0.0
LC_R2_df$sign <- LC_R2_df$V4 >= 0.0

pdf("NannoEpi_HiC_AB_log2fc.pdf", useDingbats=FALSE, width=5, height=9, family="Arial", onefile=F)

HC1_data <- HC_R1_df
LC1_data <- LC_R1_df
HC2_data <- HC_R2_df
LC2_data <- LC_R2_df

bgl_data <- bin_gene_log2fc
#str(bgl_data)

HC1_data$grp <- 'HC1'
LC1_data$grp <- 'LC1'
HC2_data$grp <- 'HC2'
LC2_data$grp <- 'LC2'

data <- rbind(HC1_data, HC2_data, LC1_data, LC2_data)
data$chr <- as.numeric(sub("chr","",data$V1))
str(data)

a <- ggdotchart(data, x="pos", y="V4", sorting="none", add = "segment", add.params=list(size=0.5,color="sign"),
           size=0.5, color="sign", dot.size=0.25, palette=c("blue","red")) + 
     facet_grid(vars(grp),vars(chr),scales="free_x",space="free_x") +
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), panel.spacing.x = unit(0, "cm"),
        strip.text.x = element_text(angle=-45,size=6), strip.background=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.y = element_text(color="black")) +
  xlab("Chromosome") + ylab("Eigenvector")

print("a")

bgl_data$HC <- bgl_data$Control
bgl_data$LC <- bgl_data$VLC
ab24 <- melt(bgl_data,measure.vars=c("HC","LC"),id.vars=c("V1","pos"),variable.name="grp")
str(ab24)
ab24_df <- merge(ab24,data,by=c("V1","pos"))
#ab24_df$value <- ab24_df$value+1
ab24_df$x <- paste(ab24_df$grp.y,ab24_df$sign,sep="_")
ab24_df$x <- factor(ab24_df$x,levels=c("HC1_TRUE","HC1_FALSE","LC1_TRUE","LC1_FALSE",
                                       "HC2_TRUE","HC2_FALSE","LC2_TRUE","LC2_FALSE"),ordered=T)
str(ab24_df)

cmps<-list(c("HC1_TRUE","HC1_FALSE"),c("LC1_TRUE","LC1_FALSE"),c("HC2_TRUE","HC2_FALSE"),c("LC2_TRUE","LC2_FALSE"))

ab <- ggbarplot(ab24_df, x="x", y="value", group="sign", fill="sign", palette=c("blue","red"), 
                add="mean_se", error.plot="errorbar", outlier.shape=NA) + coord_cartesian(ylim = c(0, 250)) +
  scale_x_discrete(labels=c("HC1_TRUE"="HC1_A","HC1_FALSE"="HC1_B","LC1_TRUE"="LC1_A","LC1_FALSE"="LC1_B",
                            "HC2_TRUE"="HC2_A","HC2_FALSE"="HC2_B","LC2_TRUE"="LC2_A","LC2_FALSE"="LC2_B")) + theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("RNA-Seq (TPM)") + 
  stat_compare_means(comparisons=cmps, method="t.test", 
                     label.y=200, tip.length=0.003)

print("ab")

data2 <- merge(HC1_data,LC1_data,by=c("V1","pos"))
df <- data2[,c("pos","sign.x","sign.y")]
df$type <- paste(data2$sign.x,data2$sign.y,sep="_")
data3 <- merge(bgl_data,df)
data3$type <- factor(data3$type, levels=c("TRUE_FALSE","FALSE_TRUE","TRUE_TRUE","FALSE_FALSE"), ordered=T)
cmps2 <- list(c("TRUE_FALSE","FALSE_TRUE"),c("TRUE_TRUE","FALSE_FALSE"))
data3$grp <- 'Rep1'

data4 <- merge(HC2_data,LC2_data,by=c("V1","pos"))
df2 <- data4[,c("pos","sign.x","sign.y")]
df2$type <- paste(data4$sign.x,data4$sign.y,sep="_")
data5 <- merge(bgl_data,df2)
data5$type <- factor(data5$type, levels=c("TRUE_FALSE","FALSE_TRUE","TRUE_TRUE","FALSE_FALSE"), ordered=T)
#cmps2 <- list(c("TRUE_FALSE","FALSE_TRUE"),c("TRUE_TRUE","FALSE_FALSE"))
data5$grp <- 'Rep2'

data35 <- rbind(data3, data5)

str(data35)

b <- ggboxplot(data35, x="type", y="V2.y", notch=T, outlier.shape=NA, bxp.errorbar=T, facet.by="grp", ncol=2) +
  scale_x_discrete(labels=c("TRUE_FALSE"="A to B","FALSE_TRUE"="B to A", "TRUE_TRUE"="A to A", "FALSE_FALSE"="B to B")) +
  theme_bw() +
  theme(text=element_text(size=11, family="Arial"), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  xlab("") + ylab("RNA-Seq (log2fc)") +
  stat_compare_means(comparisons=cmps2, method="t.test", label.y=5)

ggarrange(a,ab,b,ncol=1)

dev.off()

