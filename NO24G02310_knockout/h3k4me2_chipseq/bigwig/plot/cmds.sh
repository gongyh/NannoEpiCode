#!/bin/bash

conda activate pygenometracks

pyGenomeTracks --tracks tracks_CAH1.ini --region chr20:209500-214500 --width 20 --outFileName H3K4me2_CAH1.pdf

pyGenomeTracks --tracks tracks_ME1.ini --region chr26:168200-174300 --width 20 --outFileName H3K4me2_ME1.pdf

pyGenomeTracks --tracks tracks_TIA1.ini --region chr2:1005000-1009000 --width 20 --outFileName H3K4me2_TIA1.pdf

pyGenomeTracks --tracks tracks_MT.ini --region chr24:617500-631750 --width 30 --outFileName H3K4me2_MT.pdf

pyGenomeTracks --tracks tracks_HITs.ini --region chr14:514000-517000 --width 20 --outFileName H3K4me2_HITs.pdf

wiggletools mean WT_R1.bigWig WT_R2.bigWig WT_R3.bigWig | wigToBigWig /dev/stdin ../../../../genome/genome.fa.sizes WT_mean.bw
wiggletools mean M4_R1.bigWig M4_R2.bigWig M4_R3.bigWig | wigToBigWig /dev/stdin ../../../../genome/genome.fa.sizes M4_mean.bw
wiggletools mean M6_R1.bigWig M6_R2.bigWig M6_R3.bigWig | wigToBigWig /dev/stdin ../../../../genome/genome.fa.sizes M6_mean.bw

pyGenomeTracks --tracks tracks_HINT.ini --region chr14:514000-517000 --width 15 --outFileName H3K4me2_HINT.pdf
pyGenomeTracks --tracks tracks_PMA2.ini --region chr17:162000-168300 --width 15 --outFileName H3K4me2_PMA2.pdf

cat NO24G02310_refined.gtf | awk -F '\t' 'BEGIN{OFS="\t"}{if($4>619944){start=$4+1}else{start=$4} if($5>619944){end=$5+1}else{end=$5} print $1,$2,$3,start,end,$6,$7,$8,$9}' > NO24G02310_fix1.gtf
cat NO24G02310_fix1.gtf | awk -F '\t' 'BEGIN{OFS="\t"}{if($4>620985){start=$4+1}else{start=$4} if($5>620985){end=$5+1}else{end=$5} print $1,$2,$3,start,end,$6,$7,$8,$9}' > NO24G02310_fix2.gtf
cat NO24G02310_fix2.gtf | awk -F '\t' 'BEGIN{OFS="\t"}{if($4>627084){start=$4+1}else{start=$4} if($5>627084){end=$5+1}else{end=$5} print $1,$2,$3,start,end,$6,$7,$8,$9}' > NO24G02310_fix3.gtf
fasta_formatter -i chr24.fa -w 80 -o chr24_w80.fa
agat_sp_fix_cds_phases.pl -g NO24G02310_fix3.gtf -f chr24_w80.fa -o NO24G02310_fix4.gff
gffread -T -g chr24_w80.fa -x cds.fa -y protein.faa -o NO24G02310_fix4.gtf NO24G02310_fix4.gff 

pyGenomeTracks --tracks tracks_PMA2.ini --region chr17:160100-176300 --width 10 --height 10 --outFileName H3K4me2_PMA2_around.pdf
