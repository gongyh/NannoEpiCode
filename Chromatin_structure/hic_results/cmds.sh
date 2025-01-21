#!/bin/bash

for chr in {1..30}; do
  Rscript compare.R -i matrix/NannoH0/iced/5000/NannoH0_5000_iced.matrix -I matrix/NannoH24/iced/5000/NannoH24_5000_iced.matrix -b matrix/NannoH0/raw/5000/NannoH0_5000_abs.bed -s chr${chr}_NannoH0_5000.is30001.ids20001.insulation -S chr${chr}_NannoH24_5000.is30001.ids20001.insulation -c chr${chr} -l NannoH0 -L NannoH24 -o NannoH0_NannoH24_5000_chr${chr}.pdf -r 5000
  Rscript compare.R -i matrix/H0rep/iced/5000/H0rep_5000_iced.matrix -I matrix/H24rep/iced/5000/H24rep_5000_iced.matrix -b matrix/H0rep/raw/5000/H0rep_5000_abs.bed -s chr${chr}_H0rep_5000.is30001.ids20001.insulation -S chr${chr}_H24rep_5000.is30001.ids20001.insulation -c chr${chr} -l H0rep -L H24rep -o H0rep_H24rep_5000_chr${chr}.pdf -r 5000
done
