#!/bin/bash

export PATH=$HOME/progs/samtools-0.1.19:$HOME/progs/LACHESIS/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lib:$HOME/progs/boost_1_52_0/stage/lib:$LD_LIBRARY_PATH

## bowtie2 align, using HiC-pro
bowtie2-build IMET1v2.fasta IMET1v2

export PATH=~/bin/HiC-Pro_2.9.0/bin:~/bin/HiC-Pro_2.9.0/scripts:~/progs/HiC-Pro_2.9.0/bin/utils:$PATH
~/bin/HiC-Pro_2.9.0/bin/utils/digest_genome.py -r mboi -o IMET1_MboI.bed IMET1v2.fasta
cat IMET1v2.fasta | lengths.pl > IMET1.chrom.sizes

HiC-Pro -i rawdata -o imet1_hicpro_out3 -c config-imet1.3.txt
#HiC-Pro -i imet1_hicpro_out3/bowtie_results/bwt2 -o imet1_hicpro_out3 -c config_final.txt -s proc_hic
#cd imet1_hicpro_out3
  #~/bin/HiC-Pro_2.9.0/scripts/mapped_2hic_fragments.sh -c ../config_final.txt >> hicpro2.log
#cd ..
#HiC-Pro -i imet1_hicpro_out3/bowtie_results/bwt2 -o imet1_hicpro_out3 -c config_final.txt -s quality_checks
#cd imet1_hicpro_out3
#./plot.sh
#cd ..
#HiC-Pro -i imet1_hicpro_out3/hic_results/data -o imet1_hicpro_out3 -c config_final.txt -s merge_persample
#HiC-Pro -i imet1_hicpro_out3/hic_results/data -o imet1_hicpro_out3 -c config_final.txt -s build_contact_maps
#HiC-Pro -i imet1_hicpro_out3/hic_results/matrix -o imet1_hicpro_out3 -c config_final.txt -s ice_norm

#visualize
## Plot the genome-wide map at 1kb resolution

# 10kb
python ~/progs/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/H0rep/iced/10000/H0rep_10000_iced.matrix -o IMET1H0repgw -r 10000 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H0rep/raw/10000/H0rep_10000_ord.bed -n IMET1H0_rep -wg 1 -chr chr30 -ext pdf -dpi 300

python ~/progs/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/H24rep/iced/10000/H24rep_10000_iced.matrix -o IMET1H24repgw -r 10000 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H24rep/raw/10000/H24rep_10000_ord.bed -n IMET1H24_rep -wg 1 -chr chr30 -ext pdf -dpi 300

conda activate py27
python ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/H0rep/iced/10000/H0rep_10000_iced.matrix imet1_hicpro_out3/hic_results/matrix/H24rep/iced/10000/H24rep_10000_iced.matrix -n HC_Rep2 LC_Rep2 -wg 1 -chr chr30 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H0rep/raw/10000/H0rep_10000_ord.bed -o IMET1v2_Rep2 -r 10000 -c 1 -hR 1

python ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/10000/NannoH0_10000_iced.matrix imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/10000/NannoH24_10000_iced.matrix -n HC_Rep1 LC_Rep1 -wg 1 -chr chr30 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/NannoH0/raw/10000/NannoH0_10000_ord.bed -o IMET1v2_Rep1 -r 10000 -c 1 -hR 1

python ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/10000/NannoH0_10000_iced.matrix imet1_hicpro_out3/hic_results/matrix/H0rep/iced/10000/H0rep_10000_iced.matrix -n HC_Rep1 HC_Rep2 -wg 1 -chr chr30 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H0rep/raw/10000/H0rep_10000_ord.bed -o IMET1v2_HC -r 10000 -c 1 -hR 1

python ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/10000/NannoH24_10000_iced.matrix imet1_hicpro_out3/hic_results/matrix/H24rep/iced/10000/H24rep_10000_iced.matrix -n LC_Rep1 LC_Rep2 -wg 1 -chr chr30 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H24rep/raw/10000/H24rep_10000_ord.bed -o IMET1v2_LC -r 10000 -c 1 -hR 1

##

python ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/H0rep/iced/1000/H0rep_1000_iced.matrix imet1_hicpro_out3/hic_results/matrix/H24rep/iced/1000/H24rep_1000_iced.matrix -n HC_Rep2 LC_Rep2 -wg 0 -chr chr1 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H0rep/raw/1000/H0rep_1000_ord.bed -o IMET1v2_Rep2 -r 1000 -c 1 -hR 1

python ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/1000/NannoH0_1000_iced.matrix imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/1000/NannoH24_1000_iced.matrix -n HC_Rep1 LC_Rep1 -wg 0 -chr chr1 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/NannoH0/raw/1000/NannoH0_1000_ord.bed -o IMET1v2_Rep1 -r 1000 -c 1 -hR 1

python ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/1000/NannoH0_1000_iced.matrix imet1_hicpro_out3/hic_results/matrix/H0rep/iced/1000/H0rep_1000_iced.matrix -n HC_Rep1 HC_Rep2 -wg 0 -chr chr1 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H0rep/raw/1000/H0rep_1000_ord.bed -o IMET1v2_HC -r 1000 -c 1 -hR 1

python ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/1000/NannoH24_1000_iced.matrix imet1_hicpro_out3/hic_results/matrix/H24rep/iced/1000/H24rep_1000_iced.matrix -n LC_Rep1 LC_Rep2 -wg 0 -chr chr1 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H24rep/raw/1000/H24rep_1000_ord.bed -o IMET1v2_LC -r 1000 -c 1 -hR 1

##


python ~/src/HiCPlotter/HiCPlotter.py -chr chr1 -f imet1_hicpro_out3/hic_results/matrix/H0rep/iced/10000/H0rep_10000_iced.matrix imet1_hicpro_out3/hic_results/matrix/H24rep/iced/10000/H24rep_10000_iced.matrix -n HC LC -wg 0 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H0rep/raw/10000/H0rep_10000_ord.bed -o IMET1v2rep_chr1 -r 10000 -c 1 -hR 1 #-ptd 1 -pi 1 -fh 0

python ~/src/HiCPlotter/HiCPlotter.py -chr chr1 -f imet1_hicpro_out3/hic_results/matrix/H0rep/iced/1000/H0rep_1000_iced.matrix imet1_hicpro_out3/hic_results/matrix/H24rep/iced/1000/H24rep_1000_iced.matrix -n HC LC -wg 0 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/H0rep/raw/1000/H0rep_1000_ord.bed -o IMET1v2rep_chr1 -r 1000 -c 1 -hR 1 #-ptd 1 -pi 1 -fh 0



# use Arial font
python2 ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/1000/NannoH0_1000_iced.matrix imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/1000/NannoH24_1000_iced.matrix -n HC_Rep1 LC_Rep1 -wg 0 -chr chr1 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/NannoH0/raw/1000/NannoH0_1000_ord.bed -o IMET1v2_Rep1_arial -r 1000 -c 1 -hR 1

python2 ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/10000/NannoH0_10000_iced.matrix imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/10000/NannoH24_10000_iced.matrix -n HC_Rep1 LC_Rep1 -wg 1 -chr chr30 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/NannoH0/raw/10000/NannoH0_10000_ord.bed -o IMET1v2_Rep1_arial -r 10000 -c 1 -hR 1

python2 ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/1000/NannoH0_1000_iced.matrix imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/1000/NannoH24_1000_iced.matrix -n HC_Rep1 LC_Rep1 -wg 0 -chr chr24 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/NannoH0/raw/1000/NannoH0_1000_ord.bed -o IMET1v2_Rep1_chr24 -r 1000 -c 1 -hR 1

python2 ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/10000/NannoH0_10000_iced.matrix imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/10000/NannoH24_10000_iced.matrix -n HC_Rep1 LC_Rep1 -wg 0 -chr chr24 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/NannoH0/raw/10000/NannoH0_10000_ord.bed -o IMET1v2_Rep1_chr24 -r 10000 -c 1 -hR 1

python2 ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/1000/NannoH0_1000_iced.matrix imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/1000/NannoH24_1000_iced.matrix -n HC_Rep1 LC_Rep1 -wg 0 -chr chr14 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/NannoH0/raw/1000/NannoH0_1000_ord.bed -o IMET1v2_Rep1_chr14 -r 1000 -c 1 -hR 1

python2 ~/src/HiCPlotter/HiCPlotter.py -f imet1_hicpro_out3/hic_results/matrix/NannoH0/iced/10000/NannoH0_10000_iced.matrix imet1_hicpro_out3/hic_results/matrix/NannoH24/iced/10000/NannoH24_10000_iced.matrix -n HC_Rep1 LC_Rep1 -wg 0 -chr chr14 -ext pdf -dpi 300 -tri 1 -bed imet1_hicpro_out3/hic_results/matrix/NannoH0/raw/10000/NannoH0_10000_ord.bed -o IMET1v2_Rep1_chr14 -r 10000 -c 1 -hR 1

