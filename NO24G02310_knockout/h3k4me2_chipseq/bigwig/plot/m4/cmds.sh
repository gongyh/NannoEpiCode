#!/bin/bash

conda activate pygenometracks

pyGenomeTracks --tracks tracks_CAH1.ini --region chr20:209500-214500 --width 20 --outFileName H3K4me2_CAH1.pdf

pyGenomeTracks --tracks tracks_ME1.ini --region chr26:168200-174300 --width 20 --outFileName H3K4me2_ME1.pdf

pyGenomeTracks --tracks tracks_TIA1.ini --region chr2:1005000-1009000 --width 20 --outFileName H3K4me2_TIA1.pdf

pyGenomeTracks --tracks tracks_MT.ini --region chr24:617500-631750 --width 30 --outFileName H3K4me2_MT.pdf

cat ../NO24G02310_refined.gtf | awk -F '\t' 'BEGIN{OFS="\t"}{if($4>619900){start=$4-1}else{start=$4} if($5>619900){end=$5-1}else{end=$5} print $1,$2,$3,start,end,$6,$7,$8,$9}' > NO24G02310_fix0.gtf
cat NO24G02310_fix0.gtf | awk -F '\t' 'BEGIN{OFS="\t"}{if($4>619944){start=$4+1}else{start=$4} if($5>619944){end=$5+1}else{end=$5} print $1,$2,$3,start,end,$6,$7,$8,$9}' > NO24G02310_fix1.gtf
cat NO24G02310_fix1.gtf | awk -F '\t' 'BEGIN{OFS="\t"}{if($4>620985){start=$4+1}else{start=$4} if($5>620985){end=$5+1}else{end=$5} print $1,$2,$3,start,end,$6,$7,$8,$9}' > NO24G02310_fix2.gtf
cat NO24G02310_fix2.gtf | awk -F '\t' 'BEGIN{OFS="\t"}{if($4>627084){start=$4+1}else{start=$4} if($5>627084){end=$5+1}else{end=$5} print $1,$2,$3,start,end,$6,$7,$8,$9}' > NO24G02310_fix3.gtf
fasta_formatter -i chr24_m4.fa -w 80 -o chr24_w80.fa
agat_sp_fix_cds_phases.pl -g NO24G02310_fix3.gtf -f chr24_w80.fa -o NO24G02310_fix4.gff
gffread -T -g chr24_w80.fa -x cds.fa -y protein.faa -o NO24G02310_fix4.gtf NO24G02310_fix4.gff 
