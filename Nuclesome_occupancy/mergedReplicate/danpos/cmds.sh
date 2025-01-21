#!/bin/bash

#sortBed -i H0.mRp.clN.Fnor.smooth.positions.summit.bed > H0_summit_sorted.bed
#sortBed -i H24.mRp.clN.Fnor.smooth.positions.summit.bed > H24_summit_sorted.bed

awk -F'\t' 'BEGIN{pre_chr="";pre_pos=-1;}{if($1==pre_chr) {print $2-pre_pos} pre_chr=$1; pre_pos=$2}' H0_summit_sorted.bed > H0_phases.txt
awk -F'\t' 'BEGIN{pre_chr="";pre_pos=-1;}{if($1==pre_chr) {print $2-pre_pos} pre_chr=$1; pre_pos=$2}' H24_summit_sorted.bed > H24_phases.txt

bedtools intersect -a H0.mRp.clN.Fnor.smooth.positions.bed -b H24.mRp.clN.Fnor.smooth.positions.bed -v > HC_only.positions.bed
bedtools intersect -b H0.mRp.clN.Fnor.smooth.positions.bed -a H24.mRp.clN.Fnor.smooth.positions.bed -v > LC_only.positions.bed

export PATH=$PATH:~/src/homer/bin
perl ~/src/homer/bin/annotatePeaks.pl HC_only.positions.bed none -gtf genes.gtf > HC_only.positions.annotation

docker run -v /mnt/scc8t/gongyh/:/mnt/scc8t/gongyh/ -w $PWD nfcore/mnaseseq:dev danpos.py dpos H0/H0.mRp.clN.bed:H24/H24.mRp.clN.bed --span 1 --smooth_width 20 --width 40 --count 1000000 --out cmp --paired 1 --fdr 1

danpos.py dpos H0/H0.mRp.clN.bed:H24/H24.mRp.clN.bed --span 1 --smooth_width 20 --width 40 --count 1000000 --out cmp --paired 1 --fdr 1


bedtools makewindows -b genome.bed -w 5000 | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,"bin_"NR}' | grep chr > genome_5k.bed

bedtools intersect -a genome_5k.bed -b H0.mRp.clN.Fnor.smooth.positions.summit.bed -wao | awk -F'\t' 'BEGIN{prev="bin_0"}{bin=$4;if(bin!=prev){count[bin]=0;prev=bin} if($8!="."){count[bin]+=1}}END{for(b in count){print b"\t"count[b]}}' | sort -V > H0_summit_count_bin5k.txt

bedtools intersect -a genome_5k.bed -b H24.mRp.clN.Fnor.smooth.positions.summit.bed -wao | awk -F'\t' 'BEGIN{prev="bin_0"}{bin=$4;if(bin!=prev){count[bin]=0;prev=bin} if($8!="."){count[bin]+=1}}END{for(b in count){print b"\t"count[b]}}' | sort -V > H24_summit_count_bin5k.txt

join -1 4 -2 1 -t $'\t' genome_5k.bed H0_summit_count_bin5k.txt > H0_summit_count_bin5k.tsv
join -1 4 -2 1 -t $'\t' genome_5k.bed H24_summit_count_bin5k.txt > H24_summit_count_bin5k.tsv

