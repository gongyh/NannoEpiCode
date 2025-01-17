#!/bin/bash

python compare.py gene_5mC_HC.txt gene_5mC_LC.txt > gene_5mC_diff.txt
python compare.py genePromoter_5mC_HC.txt genePromoter_5mC_LC.txt > genePromoter_5mC_diff.txt

## Figure S1A
Rscript context2.R

## Figure S1B
Rscript drawMLgc2.R

## Figre S1C
Rscript drawDiff.R
