#!/bin/bash

#Rscript compare1k.R
#Rscript compare2k.R

cat Rep1_1k.txt Rep2_1k.txt | cut -f2-7 | grep -v IF1 | sort -k1,1V -k2,2n | uniq -d > 1k_common.txt


