#!/bin/bash

nextflow run /mnt/data7/gongyh/images/chipseq-1.2.1/main.nf -profile docker --input design2.csv --genome IMET1v2 \
  --macs_gsize 3.0e7 --save_reference --save_trimmed --save_macs_pileup --outdir results_final -resume --narrow_peak

