#!/bin/bash

res=40000

samples=(NannoH0 NannoH24 H0rep H24rep)
for sample in ${samples[*]}; do
  #contacts=../data/${sample}/${sample}_allValidPairs
  #cut -f1-7 $contacts > ${sample}.homer
  ~/src/CscoreTool/CscoreTool1.1 bins_${res}.bed ${sample}.homer ${sample}_${res} 48 $res
done

