#!/bin/bash -x
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt 
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
cat IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'
cat IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'
cat IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHH\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' 
cat IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'  
cat IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'
cat IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHH\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' 
cat IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'  
cat IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'
cat IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHH\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' 
cat IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'   
cat IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'
cat IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHH\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' 
cat IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'  
cat IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'
cat IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHH\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' 
cat IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'  
cat IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHG\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}'
cat IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | grep $'\tCHH\t' | awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>0) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' 
awk -F'\t' 'BEGIN{CG=0;CHG=0;CHH=0;OFS="\t"}{if($4>1) {if($6=="CG"){CG+=1;}else if($6=="CHG"){CHG+=1;}else if($6="CHH"){CHH+=1;}else{print $0}}}END{print CG,CHG,CHH,CG+CHG+CHH}' IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt 
awk -F'\t' 'BEGIN{CG=0;CHG=0;CHH=0;OFS="\t"}{if($4>1) {if($6=="CG"){CG+=1;}else if($6=="CHG"){CHG+=1;}else if($6="CHH"){CHH+=1;}else{print $0}}}END{print CG,CHG,CHH,CG+CHG+CHH}' IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{CG=0;CHG=0;CHH=0;OFS="\t"}{if($4>1) {if($6=="CG"){CG+=1;}else if($6=="CHG"){CHG+=1;}else if($6="CHH"){CHH+=1;}else{print $0}}}END{print CG,CHG,CHH,CG+CHG+CHH}' IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{CG=0;CHG=0;CHH=0;OFS="\t"}{if($4>1) {if($6=="CG"){CG+=1;}else if($6=="CHG"){CHG+=1;}else if($6="CHH"){CHH+=1;}else{print $0}}}END{print CG,CHG,CHH,CG+CHG+CHH}' IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{CG=0;CHG=0;CHH=0;OFS="\t"}{if($4>1) {if($6=="CG"){CG+=1;}else if($6=="CHG"){CHG+=1;}else if($6="CHH"){CHH+=1;}else{print $0}}}END{print CG,CHG,CHH,CG+CHG+CHH}' IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{CG=0;CHG=0;CHH=0;OFS="\t"}{if($4>1) {if($6=="CG"){CG+=1;}else if($6=="CHG"){CHG+=1;}else if($6="CHH"){CHH+=1;}else{print $0}}}END{print CG,CHG,CHH,CG+CHG+CHH}' IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
cat IMET1v2.*_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt | awk -F'\t' '{if($4+$5>=2) print $1,$2,$3}' | sort | uniq | wc -l
wc -l IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
