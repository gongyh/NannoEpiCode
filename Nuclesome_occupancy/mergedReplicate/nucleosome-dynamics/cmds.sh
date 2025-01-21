#!/bin/bash

#docker run -v $PWD:$PWD -v $PWD/../:$PWD/../ -w $PWD mmbirb/nucldyn readBAM --input ../H0.mRp.clN.sorted.bam --output $PWD/HC.RData.gz --type paired
#docker run -v $PWD:$PWD -v $PWD/../:$PWD/../ -w $PWD mmbirb/nucldyn readBAM --input ../H24.mRp.clN.sorted.bam --output $PWD/LC.RData.gz --type paired

#docker run -v $PWD:$PWD -v $PWD/../:$PWD/../ -w $PWD mmbirb/nucldyn nucleR --input $PWD/HC.RData.gz --output $PWD/HC_nucleR.gff --type paired
#docker run -v $PWD:$PWD -v $PWD/../:$PWD/../ -w $PWD mmbirb/nucldyn nucleR --input $PWD/LC.RData.gz --output $PWD/LC_nucleR.gff --type paired

docker run -v $PWD:$PWD -v $PWD/../:$PWD/../ -w $PWD mmbirb/nucldyn nucDyn --input1 $PWD/HC.RData.gz --input2 $PWD/LC.RData.gz --calls1 $PWD/HC_nucleR2.gff --calls2 $PWD/LC_nucleR2.gff -outputGff cmp.gff --outputBigWig cmp.bw --genome chrom.sizes --cores 32

