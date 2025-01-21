#!/bin/bash

java -mx40G -jar ~/progs/ChromHMM/ChromHMM.jar BinarizeBam -gzip -paired CHROMSIZES/IMET1v2.txt BAMS/ design.txt BinarizeBam

#mkdir -p COORDS/IMET1v2
#mkdir -p ANCHORFILES/IMET1v2

#cat IMET1v2.primary.gff3 | grep exon | cut -f1,4,5 | awk -F '\t' '{a=$2-1;print $1"\t"a"\t"$3}' | sort -k1,1 -k2,2n -k3,3n | uniq | gzip -c > COORDS/IMET1v2/RefSeqExon.IMET1v2.bed.gz

#cat IMET1v2.gene.gff3 | grep gene | cut -f1,4,5 | awk -F '\t' '{a=$2-1;print $1"\t"a"\t"$3}' | sort -k1,1 -k2,2n -k3,3n | uniq | gzip -c > COORDS/IMET1v2/RefSeqGene.IMET1v2.bed.gz

#python gff2coords.py IMET1v2.primary.gff3 CHROMSIZES/IMET1v2.txt > coords.all
#cat coords.all | grep $'TSS\t' | sort -k2,2 -k3,3n -k4,4n | uniq | cut -f2,3,4 | gzip -c > COORDS/IMET1v2/RefSeqTSS.IMET1v2.bed.gz
#cat coords.all | grep $'TES\t' | sort -k2,2 -k3,3n -k4,4n | uniq | cut -f2,3,4 | gzip -c > COORDS/IMET1v2/RefSeqTES.IMET1v2.bed.gz
#cat coords.all | grep $'TSS2kb\t' | sort -k2,2 -k3,3n -k4,4n | uniq | cut -f2,3,4 | gzip -c > COORDS/IMET1v2/RefSeqTSS2kb.IMET1v2.bed.gz
#cat coords.all | grep $'TSSanchor\t' | sort -k2,2 -k3,3n -k4,4n | uniq | cut -f2,3,4 | gzip -c > ANCHORFILES/IMET1v2/RefSeqTSS.IMET1v2.txt.gz
#cat coords.all | grep $'TESanchor\t' | sort -k2,2 -k3,3n -k4,4n | uniq | cut -f2,3,4 | gzip -c > ANCHORFILES/IMET1v2/RefSeqTES.IMET1v2.txt.gz

java -mx120G -jar ~/progs/ChromHMM/ChromHMM.jar LearnModel -p 48 -u COORDS -v ANCHORFILES -r 500 BinarizeBam Model3 3 IMET1v2
java -mx120G -jar ~/progs/ChromHMM/ChromHMM.jar LearnModel -p 48 -u COORDS -v ANCHORFILES -r 500 BinarizeBam Model4 4 IMET1v2
java -mx120G -jar ~/progs/ChromHMM/ChromHMM.jar LearnModel -p 48 -u COORDS -v ANCHORFILES -r 500 BinarizeBam Model5 5 IMET1v2
java -mx120G -jar ~/progs/ChromHMM/ChromHMM.jar LearnModel -p 48 -u COORDS -v ANCHORFILES -r 500 BinarizeBam Model6 6 IMET1v2
java -mx120G -jar ~/progs/ChromHMM/ChromHMM.jar LearnModel -p 48 -u COORDS -v ANCHORFILES -r 500 BinarizeBam Model7 7 IMET1v2
java -mx120G -jar ~/progs/ChromHMM/ChromHMM.jar LearnModel -p 48 -u COORDS -v ANCHORFILES -r 500 BinarizeBam Model8 8 IMET1v2
java -mx120G -jar ~/progs/ChromHMM/ChromHMM.jar LearnModel -p 48 -u COORDS -v ANCHORFILES -r 500 BinarizeBam Model9 9 IMET1v2

rm -rf emissions
mkdir -p emissions
cd emissions
ln -s ../Model[3-8]/emissions_*.txt .
cd ..
java -mx40G -jar ~/progs/ChromHMM/ChromHMM.jar CompareModels Model9/emissions_9.txt emissions compare


###############
python3 gff2TSS_bed.py IMET1v2.geneOnly.gff3 > IMET1v2.TSS.bed
bedtools intersect -wb -a Model5/H0_5_dense.bed -b IMET1v2.TSS.bed | cut -f4,13 | sort -k2,2 > H0_5_gene_state.txt
bedtools intersect -wb -a Model5/H24_5_dense.bed -b IMET1v2.TSS.bed | cut -f4,13 | sort -k2,2 > H24_5_gene_state.txt
python3.6 gene_state_change.py H0_5_gene_state.txt H24_5_gene_state.txt  H0_H24_CS.txt > H0_H24_CS.details


python3.6 combine.py H0_H24_CS.details HC_VLC_DEG_Up.txt HC_VLC_DEG_Down.txt > H0_H24_CS_DEG.txt


