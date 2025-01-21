#!/bin/bash

fasta=/mnt/ds28b/gongyh/Nanno/HiC/30/ordered/corrected/hicpro/IMET1v2.fasta
chrsize=/mnt/ds28b/gongyh/Nanno/HiC/30/ordered/corrected/hicpro/IMET1.chrom.sizes
res=40000
threads=48

#cooler makebins ${chrsize} ${res} > bins_${res}.bed
#cooltools genome binnify --all-names ${chrsize} ${res} > genome_bins_${res}.txt
#cooltools genome gc genome_bins_${res}.txt ${fasta} > genome_gc_${res}.txt

samples=(NannoH0 NannoH24 H0rep H24rep)
for sample in ${samples[*]}; do
  #contacts=../data/${sample}/${sample}_allValidPairs
  #cooler cload pairs -c1 2 -p1 3 -c2 5 -p2 6 bins_${res}.bed ${contacts} ${sample}_${res}.cool
  echo "${sample}: step 1 done"
  #cp ${sample}_${res}.cool ${sample}_${res}_norm.cool
  #cooler balance ${sample}_${res}_norm.cool -p ${threads} --force
  echo "${sample}: step 2 done"
  #cooltools eigs-cis --phasing-track bins_${res}_gd.txt --out-prefix ${sample}_${res} ${sample}_${res}_norm.cool
  #cut -f1,2,3,5 ${sample}_${res}.cis.vecs.tsv | grep -v E1 > ${sample}_${res}.cis.E1.bedgraph
done

## HC
cat ../data/NannoH0/NannoH0_allValidPairs ../data/H0rep/H0rep_allValidPairs > HC_allValidPairs
cooler cload pairs -c1 2 -p1 3 -c2 5 -p2 6 bins_${res}.bed HC_allValidPairs HC_${res}.cool
cp HC_${res}.cool HC_${res}_norm.cool
cooler balance HC_${res}_norm.cool -p ${threads} --force
cooltools eigs-cis --phasing-track bins_${res}_gd.txt --out-prefix HC_${res} HC_${res}_norm.cool
cut -f1,2,3,5 HC_${res}.cis.vecs.tsv | grep -v E1 > HC_${res}.cis.E1.bedgraph

cooltools expected-cis -c ${res} -o HC_exp.tsv HC_${res}_norm.cool
cooltools saddle --contact-type cis --strength --qrange 0.025 0.975 -o chr1_saddle --fig pdf HC_${res}_norm.cool HC_${res}.cis.vecs.tsv HC_exp.tsv

## LC
cat ../data/NannoH24/NannoH24_allValidPairs ../data/H24rep/H24rep_allValidPairs > LC_allValidPairs
cooler cload pairs -c1 2 -p1 3 -c2 5 -p2 6 bins_${res}.bed LC_allValidPairs LC_${res}.cool
cp LC_${res}.cool LC_${res}_norm.cool
cooler balance LC_${res}_norm.cool -p ${threads} --force
cooltools eigs-cis --phasing-track bins_${res}_gd.txt --out-prefix LC_${res} LC_${res}_norm.cool
cut -f1,2,3,5 LC_${res}.cis.vecs.tsv | grep -v E1 > LC_${res}.cis.E1.bedgraph

cooltools expected-cis --view chr1.bed -c ${res} -o LC_exp.tsv LC_${res}_norm.cool

