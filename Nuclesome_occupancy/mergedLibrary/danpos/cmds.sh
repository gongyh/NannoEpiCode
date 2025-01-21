#!/bin/bash

docker run -v /mnt/scc8t/gongyh/:/mnt/scc8t/gongyh/ -w $PWD nfcore/mnaseseq:dev danpos.py dpos H0/:H24/ \
  --span 1 --smooth_width 20 --width 40 --count 1000000 --out cmp --paired 1 --fdr 1

danpos.py dpos H0/:H24/ --span 1 --smooth_width 20 --width 40 --count 1000000 --out cmp --paired 1 --fdr 1

cat H*.Fnor.smooth.positions.bed | sort -k1,1 -k2,2n | bedtools merge -i - > consensus_peaks.bed
echo $'GeneID\tChr\tStart\tEnd\tStrand' > consensus_peaks.saf
awk '{print $1":"$2"-"$3"\t"$1"\t"$2"\t"$3"\t.\t+"}' consensus_peaks.bed >> consensus_peaks.saf
#bedtools multicov -bams ../H0_R1.mLb.clN.sorted.bam  ../H0_R2.mLb.clN.sorted.bam  ../H24_R1.mLb.clN.sorted.bam  ../H24_R2.mLb.clN.sorted.bam -bed consensus_peaks.bed > peak_counts.txt
featureCounts -a consensus_peaks.saf -o peak_counts.txt -F SAF -p -O --fracOverlap 0.2 --donotsort -T 48 ../H0_R1.mLb.clN.sorted.bam  ../H0_R2.mLb.clN.sorted.bam  ../H24_R1.mLb.clN.sorted.bam  ../H24_R2.mLb.clN.sorted.bam

cat peak_counts.txt | awk -F'\t' 'BEGIN{OFS="\t"}NR>2{print $1,$2,$3,$4,"+"}' > mnase.consense_peaks.txt
annotatePeaks.pl mnase.consense_peaks.txt genome.fa -gid -gtf genes.gtf -cpu 6  > mnase.consense_peaks.annotatePeaks.txt

