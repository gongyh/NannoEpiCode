#!/bin/bash

TMP=$PWD/tmp
TMPDIR=$TMP
TEMPDIR=$TMP
export TMP TMPDIR TEMPDIR

export _JAVA_OPTIONS="-Djava.io.tmpdir=$PWD/tmp"

#~/bin/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/NannoH0/H0_Lib1_Lane1_IMET1v2.bwt2pairs.validPairs -g ../IMET1.chrom.sizes -j ~/progs/juicer-1.5.6/juicer_tools_linux_0.8.jar -r ../IMET1_MboI.bed -t tmp -o juicer

####java -jar /home/gene/gongyh/progs/juicer-1.5.6/juicer_tools_linux_0.8.jar pre -f tmp/39017_resfrag.juicebox tmp/39017_allValidPairs.pre_juicebox_sorted juicer/H0Rep1.hic ../IMET1.chrom.sizes

~/bin/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/NannoH24/H24_Lib1_Lane12_IMET1v2.bwt2pairs.validPairs -g ../IMET1.chrom.sizes -j ~/progs/juicer-1.5.6/juicer_tools_linux_0.8.jar -r ../IMET1_MboI.bed -t tmp -o juicer

~/bin/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/H0rep/0h_100_IMET1v2.bwt2pairs.validPairs -g ../IMET1.chrom.sizes -j ~/progs/juicer-1.5.6/juicer_tools_linux_0.8.jar -r ../IMET1_MboI.bed -t tmp -o juicer

~/bin/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/H24rep/24h_100_IMET1v2.bwt2pairs.validPairs -g ../IMET1.chrom.sizes -j ~/progs/juicer-1.5.6/juicer_tools_linux_0.8.jar -r ../IMET1_MboI.bed -t tmp -o juicer

