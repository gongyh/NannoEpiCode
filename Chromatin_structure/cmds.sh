#!/bin/bash

# diffHic to filter artifactual read pairs
# mark duplicates using Picard
#java -jar ~/progs/picard-tools-1.139/picard.jar SortSam I=bowtie_results/bwt2/NannoH0/H0_Lib1_Lane1_IMET1v2.bwt2pairs.bam O=NannoH0_sort.bam SORT_ORDER=coordinate TMP_DIR=tmp
#java -jar ~/progs/picard-tools-1.139/picard.jar MarkDuplicates I=NannoH0_sort.bam O=NannoH0_md.bam M=NannoH0_md.txt TMP_DIR=tmp
#java -Xmx20G -jar ~/progs/picard-tools-1.139/picard.jar SortSam I=NannoH0_md.bam O=NannoH0_md_sort.bam SORT_ORDER=queryname TMP_DIR=tmp2

#java -jar ~/progs/picard-tools-1.139/picard.jar SortSam I=bowtie_results/bwt2/NannoH24/H24_Lib1_Lane12_IMET1v2.bwt2pairs.bam O=NannoH24_sort.bam SORT_ORDER=coordinate TMP_DIR=tmp
#java -jar ~/progs/picard-tools-1.139/picard.jar MarkDuplicates I=NannoH24_sort.bam O=NannoH24_md.bam M=NannoH24_md.txt TMP_DIR=tmp
#java -Xmx200G -jar ~/progs/picard-tools-1.139/picard.jar SortSam I=NannoH24_md.bam O=NannoH24_md_sort.bam SORT_ORDER=queryname TMP_DIR=tmp

'''
java -Xmx200G -jar ~/progs/picard-tools-1.139/picard.jar SortSam I=bowtie_results/bwt2/H0rep/0h_100_IMET1v2.bwt2pairs.bam O=NannoH0rep_sort.bam SORT_ORDER=coordinate TMP_DIR=tmp
java -Xmx200G -jar ~/progs/picard-tools-1.139/picard.jar MarkDuplicates I=NannoH0rep_sort.bam O=NannoH0rep_md.bam M=NannoH0rep_md.txt TMP_DIR=tmp
java -Xmx200G -jar ~/progs/picard-tools-1.139/picard.jar SortSam I=NannoH0rep_md.bam O=NannoH0rep_md_sort.bam SORT_ORDER=queryname TMP_DIR=tmp

java -Xmx200G -jar ~/progs/picard-tools-1.139/picard.jar SortSam I=bowtie_results/bwt2/H24rep/24h_100_IMET1v2.bwt2pairs.bam O=NannoH24rep_sort.bam SORT_ORDER=coordinate TMP_DIR=tmp
java -Xmx200G -jar ~/progs/picard-tools-1.139/picard.jar MarkDuplicates I=NannoH24rep_sort.bam O=NannoH24rep_md.bam M=NannoH24rep_md.txt TMP_DIR=tmp
java -Xmx200G -jar ~/progs/picard-tools-1.139/picard.jar SortSam I=NannoH24rep_md.bam O=NannoH24rep_md_sort.bam SORT_ORDER=queryname TMP_DIR=tmp
'''

'''
python filterDupAnd4k.py NannoH0_md_sort.bam NannoH0_filter.bam
python filterDupAnd4k.py NannoH24_md_sort.bam NannoH24_filter.bam
python filterDupAnd4k.py NannoH0rep_md_sort.bam NannoH0rep_filter.bam
python filterDupAnd4k.py NannoH24rep_md_sort.bam NannoH24rep_filter.bam
'''

## convert hic-pro to TADbit .tsv
# need mapped length of every reads, fragment file, and allvalid reads
## H0rep
#samtools view -f 64 bowtie_results/bwt2/H0rep/0h_100_IMET1v2.bwt2pairs.bam | awk -F '\t' '{print $1"\t"length($10)}' > H0rep_r1_map_len.txt
#samtools view -f 128 bowtie_results/bwt2/H0rep/0h_100_IMET1v2.bwt2pairs.bam | awk -F '\t' '{print $1"\t"length($10)}' > H0rep_r2_map_len.txt
## H24rep
#samtools view -f 64 bowtie_results/bwt2/H24rep/24h_100_IMET1v2.bwt2pairs.bam | awk -F '\t' '{print $1"\t"length($10)}' > H24rep_r1_map_len.txt
#samtools view -f 128 bowtie_results/bwt2/H24rep/24h_100_IMET1v2.bwt2pairs.bam | awk -F '\t' '{print $1"\t"length($10)}' > H24rep_r2_map_len.txt
## NannoH0
#samtools view -f 64 bowtie_results/bwt2/NannoH0/H0_Lib1_Lane1_IMET1v2.bwt2pairs.bam | awk -F '\t' '{print $1"\t"length($10)}' > H0_r1_map_len.txt
#samtools view -f 128 bowtie_results/bwt2/NannoH0/H0_Lib1_Lane1_IMET1v2.bwt2pairs.bam | awk -F '\t' '{print $1"\t"length($10)}' > H0_r2_map_len.txt
## NannoH24
#samtools view -f 64 bowtie_results/bwt2/NannoH24/H24_Lib1_Lane12_IMET1v2.bwt2pairs.bam | awk -F '\t' '{print $1"\t"length($10)}' > H24_r1_map_len.txt
#samtools view -f 128 bowtie_results/bwt2/NannoH24/H24_Lib1_Lane12_IMET1v2.bwt2pairs.bam | awk -F '\t' '{print $1"\t"length($10)}' > H24_r2_map_len.txt

## covert to tsv format
#python hicpro2tadbit.py hic_results/data/H0rep/H0rep_allValidPairs ../IMET1_MboI.bed H0rep_r1_map_len.txt H0rep_r2_map_len.txt H0rep.tsv
#python hicpro2tadbit.py hic_results/data/H24rep/H24rep_allValidPairs ../IMET1_MboI.bed H24rep_r1_map_len.txt H24rep_r2_map_len.txt H24rep.tsv
#python hicpro2tadbit.py hic_results/data/NannoH0/NannoH0_allValidPairs ../IMET1_MboI.bed H0_r1_map_len.txt H0_r2_map_len.txt H0.tsv
#python hicpro2tadbit.py hic_results/data/NannoH24/NannoH24_allValidPairs ../IMET1_MboI.bed H24_r1_map_len.txt H24_r2_map_len.txt H24.tsv

python drawCP.py -b hic_results/matrix/NannoH0/raw/1000/NannoH0_1000_abs.bed -o NannoH0_1000_cp.pdf hic_results/matrix/NannoH0/iced/1000/NannoH0_1000_iced.matrix
python drawCP.py -b hic_results/matrix/NannoH24/raw/1000/NannoH24_1000_abs.bed -o NannoH24_1000_cp.pdf hic_results/matrix/NannoH24/iced/1000/NannoH24_1000_iced.matrix
python drawCP.py -b hic_results/matrix/H0rep/raw/1000/H0rep_1000_abs.bed -o H0rep_1000_cp.pdf hic_results/matrix/H0rep/iced/1000/H0rep_1000_iced.matrix
python drawCP.py -b hic_results/matrix/H24rep/raw/1000/H24rep_1000_abs.bed -o H24rep_1000_cp.pdf hic_results/matrix/H24rep/iced/1000/H24rep_1000_iced.matrix

hicConvertFormat --matrices hic_results/matrix/NannoH0/iced/rfbin/NannoH0_rfbin_iced.matrix --outFileName NannoH0_chr18.cool --inputFormat hicpro --outputFormat cool --chromosome chr18 --bedFileHicpro hic_results/matrix/NannoH0/raw/rfbin/NannoH0_rfbin_abs.bed

hicConvertFormat --matrices hic_results/matrix/NannoH24/iced/rfbin/NannoH24_rfbin_iced.matrix --outFileName NannoH24_chr18.cool --inputFormat hicpro --outputFormat cool --chromosome chr18 --bedFileHicpro hic_results/matrix/NannoH24/raw/rfbin/NannoH24_rfbin_abs.bed

python2 ~/src/HiCPlotter/HiCPlotter.py -f hic_results/matrix/NannoH0/iced/40000/NannoH0_40000_iced.matrix -n HC_Rep1 -wg 0 -chr chr1 -ext pdf -dpi 300 -tri 1 -bed hic_results/matrix/NannoH0/raw/40000/NannoH0_40000_ord.bed -o IMET1v2_Rep1_chr1 -r 40000 -hR 1
python2 ~/src/HiCPlotter/HiCPlotter.py -f hic_results/matrix/NannoH0/raw/40000/NannoH0_40000.matrix -n HC_Rep1 -wg 0 -chr chr1 -ext pdf -dpi 300 -tri 1 -bed hic_results/matrix/NannoH0/raw/40000/NannoH0_40000_ord.bed -o IMET1v2_Rep1 -r 40000 -hR 1

hicConvertFormat -m hic_results/matrix20220518/HC_40000_norm.cool -o HC_40000.h5 --inputFormat cool --outputFormat h5 -r 40000
hicTransform -m HC_40000.h5 --method obs_exp -o HC_40000_obs_exp.h5
hicTransform -m HC_40000_obs_exp.h5 --method pearson -o HC_40000_pearson.h5
#hicPCA -m HC_40000.h5 -o pca1.bw pca2.bw --format bigwig
hicPlotMatrix -m HC_40000_obs_exp.h5 --outFileName HC_40000_obsexp_chr1.pdf --region chr1:0-1690242 --bigwig HC_40000.cis.E1.bw
hicPlotMatrix -m HC_40000_pearson.h5 --outFileName HC_40000_pearson_chr1.pdf --region chr1:0-1690242 --bigwig HC_40000.cis.E1.bw


