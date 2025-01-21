#!/bin/bash

samples=("NannoH0" "NannoH24" "H0rep" "H24rep")

# 1000, 2000, rfbin

res=1000

for sample in ${samples[*]}; do
  mkdir -p ${sample}_${res}
  python ~/progs/fithic-v.1.1.3/utils/hicpro2fithic.py -i ../hic_results/matrix/$sample/raw/$res/${sample}_${res}.matrix \
    -b ../hic_results/matrix/$sample/raw/$res/${sample}_${res}_ord.bed \
    -s ../hic_results/matrix/$sample/iced/$res/${sample}_${res}_iced.matrix.biases \
    -r $res -o ${sample}_${res}

done


## fithic -f fithic.fragmentMappability.gz -i fithic.interactionCounts.gz -t fithic.biases.gz -v -l H0Rep1

cat H0Rep12.short | grep -v fragmentMid1 | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$2+1;print $3,$4,$4+1}' > H0Rep12.bed
cat H24Rep12.short | grep -v fragmentMid1 | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$2+1;print $3,$4,$4+1}' > H24Rep12.bed 

bedtools intersect -a H0Rep12.bed -b IMET1v2.geneOnly.bed -wo | cut -f1,2,3,7 > H0Rep12_genes.bed
bedtools intersect -a H24Rep12.bed -b IMET1v2.geneOnly.bed -wo | cut -f1,2,3,7 > H24Rep12_genes.bed

cat H0Rep12_genes.bed H24Rep12_genes.bed | cut -f4 | sort | uniq > ovl_genes.txt

################
awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="+"){if($2<1000){newS=0}else{newS=$2-1000} print $1,newS,$3,$4}else{print $1,$2,$3+1000,$4}}' IMET1v2.geneOnly.bed > IMET1v2_gp1k.bed
cat H0Rep12.bed H24Rep12.bed > t1.bed
bedtools intersect -a rfbin_ord.bed -b t1.bed -wa | sort | uniq > t2.bed

bedtools intersect -a IMET1v2_gp1k.bed -b t2.bed -wa | sort | uniq | cut -f4 > t3.gids
grep -f t3.gids H0_H24_DEG.txt | grep -v Not | cut -f1 > g11

awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="+"){if($2<1000){newS=0}else{newS=$2-1000} print $1,newS,$2,$4}else{print $1,$3,$3+1000,$4}}' IMET1v2.geneOnly.bed > IMET1v2_promoter.bed
bedtools intersect -a IMET1v2_promoter.bed -b t2.bed -wo | sort | uniq > g10

pyGenomeTracks --tracks contacts.ini --region chr1:1674500-1681000 --width 10 --height 4 --outFileName c1.pdf
pyGenomeTracks --tracks contacts.ini --region chr6:1312000-1345000 --width 10 --height 4 --outFileName c2.pdf
pyGenomeTracks --tracks contacts.ini --region chr7:1-8000 --width 10 --height 4 --outFileName c3.pdf
pyGenomeTracks --tracks contacts.ini --region chr7:1318000-1332000 --width 10 --height 4 --outFileName c4.pdf
pyGenomeTracks --tracks contacts.ini --region chr10:1000-7000 --width 10 --height 4 --outFileName c5.pdf
pyGenomeTracks --tracks contacts.ini --region chr21:829000-835000 --width 10 --height 4 --outFileName c6.pdf

