#!/bin/bash

nextflow run /mnt/data7/gongyh/images/mnaseseq/main.nf -profile docker --input design.csv --genome IMET1v2 \
  --save_reference --save_trimmed --outdir results -resume

