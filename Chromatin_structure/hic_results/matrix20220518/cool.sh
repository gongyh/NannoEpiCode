#!/bin/bash

fasta=/mnt/ds28b/gongyh/Nanno/HiC/30/ordered/corrected/hicpro/IMET1v2.fasta
chrsize=/mnt/ds28b/gongyh/Nanno/HiC/30/ordered/corrected/hicpro/IMET1.chrom.sizes
res=10000
#sample=NannoH0
threads=48

cooler makebins ${chrsize} ${res} > bins_${res}.bed
cooltools genome binnify --all-names ${chrsize} ${res} > genome_bins_${res}.txt
cooltools genome gc genome_bins_${res}.txt ${fasta} > genome_gc_${res}.txt

samples=(NannoH24 H0rep H24rep)
samples=(NannoH0)
for sample in ${samples[*]}; do
  contacts=../data/${sample}/${sample}_allValidPairs
  cooler cload pairs -c1 2 -p1 3 -c2 5 -p2 6 bins_${res}.bed ${contacts} ${sample}_${res}.cool
  cp ${sample}_${res}.cool ${sample}_${res}_norm.cool
  cooler balance ${sample}_${res}_norm.cool -p ${threads} --force
  cooltools call-compartments --contact-type cis -o ${sample}_compartments ${sample}_${res}_norm.cool
  awk -F"\t" 'NR>1{OFS="\t"; if($4==""){$4=0}; print $1,$2,$3,$4}' ${sample}_compartments.cis.vecs.tsv | sort -k1,1 -k2,2n > ${sample}_compartments.cis.E1.bedgraph
done

