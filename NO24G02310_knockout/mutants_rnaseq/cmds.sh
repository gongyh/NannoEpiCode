#!/bin/bash

/opt/images/bin/nextflow run /mnt/data4/gongyh/Nanno/rnaseq-1.4.2/main.nf \
  --reads '../*_[12].clean.fq.gz' -profile docker \
  --genome IMET1v2 --saveReference --skipBiotypeQC \
  --aligner star --pseudo_aligner salmon \
  --outdir results -w work -resume

export TRINITY_HOME=$HOME/src/trinityrnaseq-Trinity-v2.4.0

ls results/stringtieFPKM/*gene_abund.txt > genes.listing_target_files.txt

perl $TRINITY_HOME/util/abundance_estimates_to_matrix2.pl --est_method rsem2 \
  --quant_files genes.listing_target_files.txt --out_prefix HM_genes

head -n 1 HM_genes.TMM.EXPR.matrix | sed 's/_1.cleanAligned.sortedByCoord.out.gene_abund.txt//g' > HM_genes.TMM.TPM.matrix
tail -n+2 HM_genes.TMM.EXPR.matrix | sort -k1,1 >> HM_genes.TMM.TPM.matrix

## DE analysis

$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
  --matrix results/featureCounts/merged_gene_counts.txt --method edgeR \
  --samples_file design.txt --contrasts constrasts.txt --output edgeR

cd edgeR
#$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
#  --matrix ../HM_genes.TMM.TPM.matrix -P 0.001 -C 2 --samples ../design2.txt
$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
  --matrix ../HM_genes.TMM.TPM.matrix -P 0.001 -C 1 --samples ../design2.txt
cd ..


## enrichment analysis
cut -f1 edgeR/merged_gene_counts.txt.M4_vs_WT.edgeR.DE_results.P0.001_C1.DE.subset | grep -v sampleA > M4_DEs.txt
cut -f1 edgeR/merged_gene_counts.txt.M6_vs_WT.edgeR.DE_results.P0.001_C1.DE.subset | grep -v sampleA > M6_DEs.txt
cat edgeR/merged_gene_counts.txt.M4_vs_WT.edgeR.DE_results | cut -f1,4 | grep -v sampleA > M4_logfc.tsv
cat edgeR/merged_gene_counts.txt.M6_vs_WT.edgeR.DE_results | cut -f1,4 | grep -v sampleA > M6_logfc.tsv
Rscript go_clusterProfiler4.R
Rscript go_clusterProfiler6.R

