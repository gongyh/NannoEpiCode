#!/bin/bash

awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>3) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>3) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>3) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>3) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt 
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>3) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
awk -F'\t' 'BEGIN{line=0;num=0;}{if($4+$5>3) {num+=$4/($4+$5);line+=1;}}END{print 100*num/line}' IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
