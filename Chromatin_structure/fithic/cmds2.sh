#!/bin/bash

samples=("NannoH0" "NannoH24" "H0rep" "H24rep")
res=rfbin

for sample in ${samples[*]}; do
  mkdir -p ${sample}_${res}
  python ~/progs/fithic-v.1.1.3/utils/hicpro2fithic.py -i ../hic_results/matrix/$sample/raw/$res/${sample}_${res}.matrix \
    -b ../hic_results/matrix/$sample/raw/$res/${sample}_${res}_ord.bed \
    -s ../hic_results/matrix/$sample/iced/$res/${sample}_${res}_iced.matrix.biases \
    -r 0 -o ${sample}_${res}

done


## fithic -f fithic.fragmentMappability.gz -i fithic.interactionCounts.gz -t fithic.biases.gz -o H0Rep1 -v -l H0Rep1
cd NannoH0_rfbin
gzip -cd H0Rep1/H0Rep1.spline_pass2.significances.txt.gz | awk -F '\t' 'NR==1{print $0}NR>1{if($7<=0.05) print $0}' > H0Rep1.spline_pass2.significances.FDR5.txt
cd ..
cd NannoH24_rfbin
gzip -cd H24Rep1/H24Rep1.spline_pass2.significances.txt.gz | awk -F '\t' 'NR==1{print $0}NR>1{if($7<=0.05) print $0}' > H24Rep1.spline_pass2.significances.FDR5.txt
cd ..
cd H0rep_rfbin
gzip -cd H0Rep2/H0Rep2.spline_pass2.significances.txt.gz | awk -F '\t' 'NR==1{print $0}NR>1{if($7<=0.05) print $0}' > H0Rep2.spline_pass2.significances.FDR5.txt
cd ..
cd H24rep_rfbin
gzip -cd H24Rep2/H24Rep2.spline_pass2.significances.txt.gz | awk -F '\t' 'NR==1{print $0}NR>1{if($7<=0.05) print $0}' > H24Rep2.spline_pass2.significances.FDR5.txt
cd ..

