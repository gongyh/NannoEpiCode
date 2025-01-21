#!/usr/bin/env python3.6

import sys
import math

if len(sys.argv)<3:
    print("Usage: python3.6 combineRep.py H0_H3K4me2_R1_summits_annotations.txt H0_H3K4me2_R2_summits_annotations.txt")

Rep1 = sys.argv[1]
Rep2 = sys.argv[2]

header = True
Rep1_dict = dict()
with open(Rep1) as fh:
    for line in fh:
        if header:
            header = False
            continue
        cline = line.strip().split("\t")
        if cline[13] in Rep1_dict:
            Rep1_dict[cline[13]] = min(Rep1_dict[cline[13]], int(cline[14]))
        else:
            Rep1_dict[cline[13]] = int(cline[14])

header = True
Rep2_dict = dict()
with open(Rep2) as fh:
    for line in fh:
        if header:
            header = False
            continue
        cline = line.strip().split("\t")
        if cline[13] in Rep2_dict:
            Rep2_dict[cline[13]] = min(Rep2_dict[cline[13]], int(cline[14]))
        else:
            Rep2_dict[cline[13]] = int(cline[14])

genes_list = list(Rep1_dict)
genes_list.extend(list(Rep2_dict))
genes = set(genes_list)

merged_dict = dict()
for gene in genes:
    if gene in Rep1_dict.keys():
        if gene in Rep2_dict.keys(): # average
            merged_dict[gene] = math.ceil((Rep1_dict[gene] + Rep2_dict[gene])/2)
        else: # only in Rep1
            merged_dict[gene] = Rep1_dict[gene]
    elif gene in Rep2_dict.keys(): # only in Rep2
        merged_dict[gene] = Rep2_dict[gene]
    print(gene+"\t%d"%(merged_dict[gene]))
