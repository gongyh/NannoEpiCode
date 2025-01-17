#!/bin/bash

samples=("0_1" "0_2" "0_3" "24_1" "24_2" "24_3")

for sample in ${samples[*]}; do

samtools sort IMET1v2.${sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam -@ 4 -o IMET1v2.${sample}_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bam
~/progs/qualimap_v2.1.3/qualimap bamqc -bam IMET1v2.${sample}_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bam -outdir IMET1v2.${sample} --collect-overlap-pairs -nt 48

samtools sort IMET1v2.${sample}_R1_val_1_bismark_bt2_pe.bam -@ 4 -o IMET1v2.${sample}_R1_val_1_bismark_bt2_pe.sorted.bam
~/progs/preseq_v2.0/preseq lc_extrap -v -B IMET1v2.${sample}_R1_val_1_bismark_bt2_pe.sorted.bam -o IMET1v2.${sample}.lc.txt

done


python combineCXreport.py IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > IMET1v2.H_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt

python combineCXreport.py IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > IMET1v2.L_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt

python calcMLdistr.py IMET1v2.H_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt LH_genes.TMM.TPM.avr.matrix ../IMET1v2_bt2/genes.primary.bed H ../IMET1v2_bt2/IMET1.bed > Hdistr.txt

python calcMLdistr.py IMET1v2.L_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt LH_genes.TMM.TPM.avr.matrix ../IMET1v2_bt2/genes.primary.bed L ../IMET1v2_bt2/IMET1.bed > Ldistr.txt

python calcMLdistr_fix.py IMET1v2.H_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt LH_genes.TMM.TPM.avr.matrix ../IMET1v2_bt2/genes.primary.bed H ../IMET1v2_bt2/IMET1.bed > Hdistr_fix.txt

python calcMLdistr_fix.py IMET1v2.L_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt LH_genes.TMM.TPM.avr.matrix ../IMET1v2_bt2/genes.primary.bed L ../IMET1v2_bt2/IMET1.bed > Ldistr_fix.txt

python calcMLdistr_fix2.py IMET1v2.H_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt LH_genes.TMM.TPM.avr.matrix ../IMET1v2_bt2/genes.primary.bed H ../IMET1v2_bt2/IMET1.bed > Hdistr_fix2.txt

python calcMLdistr_fix2.py IMET1v2.L_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt LH_genes.TMM.TPM.avr.matrix ../IMET1v2_bt2/genes.primary.bed L ../IMET1v2_bt2/IMET1.bed > Ldistr_fix2.txt

python3 calcMLgc_each.py ../IMET1v2_bt2/genes.primary.bed IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > gene_5mC_HC.txt

python3 calcMLgc_each.py ../IMET1v2_bt2/genes.primary.bed IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > gene_5mC_LC.txt

python3 calcMLgc_promoter.py ../IMET1v2_bt2/genes.primary.bed IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > genePromoter_5mC_HC.txt

python3 calcMLgc_promoter.py ../IMET1v2_bt2/genes.primary.bed IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > genePromoter_5mC_LC.txt

python calcMLgc_10k.py 10kb.bed IMET1v2.0_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > HC_R1_5mC_10k.txt
python calcMLgc_10k.py 10kb.bed IMET1v2.0_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > HC_R2_5mC_10k.txt
python calcMLgc_10k.py 10kb.bed IMET1v2.0_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > HC_R3_5mC_10k.txt
python calcMLgc_10k.py 10kb.bed IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > LC_R1_5mC_10k.txt
python calcMLgc_10k.py 10kb.bed IMET1v2.24_2_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > LC_R2_5mC_10k.txt
python calcMLgc_10k.py 10kb.bed IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt > LC_R3_5mC_10k.txt
