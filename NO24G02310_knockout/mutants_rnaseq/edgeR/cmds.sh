#!/bin/bash

cut -f1 merged_gene_counts.txt.M4_vs_WT.edgeR.DE_results.P0.001_C1.WT-UP.subset | sed 's/sampleA/DEG_Down/g' > M4_DEG_Down.txt
cut -f1 merged_gene_counts.txt.M4_vs_WT.edgeR.DE_results.P0.001_C1.M4-UP.subset  | sed 's/sampleA/DEG_Up/g' > M4_DEG_Up.txt
cut -f1 merged_gene_counts.txt.M6_vs_WT.edgeR.DE_results.P0.001_C1.WT-UP.subset | sed 's/sampleA/DEG_Down/g' > M6_DEG_Down.txt
cut -f1 merged_gene_counts.txt.M6_vs_WT.edgeR.DE_results.P0.001_C1.M6-UP.subset  | sed 's/sampleA/DEG_Up/g' > M6_DEG_Up.txt

cut -f1 merged_gene_counts.txt.M4_vs_WT.edgeR.DE_results.P0.001_C0.5.WT-UP.subset | sed 's/sampleA/DEG_Down/g' > M4_DEG_Down2.txt
cut -f1 merged_gene_counts.txt.M4_vs_WT.edgeR.DE_results.P0.001_C0.5.M4-UP.subset  | sed 's/sampleA/DEG_Up/g' > M4_DEG_Up2.txt
cut -f1 merged_gene_counts.txt.M6_vs_WT.edgeR.DE_results.P0.001_C0.5.WT-UP.subset | sed 's/sampleA/DEG_Down/g' > M6_DEG_Down2.txt
cut -f1 merged_gene_counts.txt.M6_vs_WT.edgeR.DE_results.P0.001_C0.5.M6-UP.subset  | sed 's/sampleA/DEG_Up/g' > M6_DEG_Up2.txt

awk -F'\t' 'BEGIN{print "DEG_Up"}NR>1{if($4>0) print $1}' merged_gene_counts.txt.M4_vs_WT.edgeR.DE_results > M4_DEG_Up3.txt
awk -F'\t' 'BEGIN{print "DEG_Down"}NR>1{if($4<0) print $1}' merged_gene_counts.txt.M4_vs_WT.edgeR.DE_results > M4_DEG_Down3.txt
awk -F'\t' 'BEGIN{print "DEG_Up"}NR>1{if($4>0) print $1}' merged_gene_counts.txt.M6_vs_WT.edgeR.DE_results > M6_DEG_Up3.txt
awk -F'\t' 'BEGIN{print "DEG_Down"}NR>1{if($4<0) print $1}' merged_gene_counts.txt.M6_vs_WT.edgeR.DE_results > M6_DEG_Down3.txt

cut -f1,4 merged_gene_counts.txt.M4_vs_WT.edgeR.DE_results | sed 's/sampleA\tlogCPM/Gene\tlogFC/g' > M4vsWT_logFC.txt
cut -f1,4 merged_gene_counts.txt.M6_vs_WT.edgeR.DE_results | sed 's/sampleA\tlogCPM/Gene\tlogFC/g' > M6vsWT_logFC.txt

csvjoin -c gene M4_vs_WT.edgeR.P0.001_C1.txt IMET1_ahrd.tsv > M4_vs_WT.edgeR.P0.001_C1.csv
csvjoin -c gene M6_vs_WT.edgeR.P0.001_C1.txt IMET1_ahrd.tsv > M6_vs_WT.edgeR.P0.001_C1.csv 

cat ../M4_DEs.txt ../M6_DEs.txt  | sort | uniq -d > M4M6_DE_common176.gids
head -n 1 M4_vs_WT.edgeR.P0.001_C1.csv > M4_vs_WT.edgeR.P0.001_C1_common.csv
head -n 1 M6_vs_WT.edgeR.P0.001_C1.csv > M6_vs_WT.edgeR.P0.001_C1_common.csv
grep -f M4M6_DE_common176.gids M4_vs_WT.edgeR.P0.001_C1.csv >> M4_vs_WT.edgeR.P0.001_C1_common.csv
grep -f M4M6_DE_common176.gids M6_vs_WT.edgeR.P0.001_C1.csv >> M6_vs_WT.edgeR.P0.001_C1_common.csv

