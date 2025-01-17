#!/usr/bin/env python

## calc methylation level of a class of genomic component

from __future__ import print_function
import sys
from multiprocessing import Pool

if len(sys.argv) != 3:
    print("Usage: python calcMLgc.py IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt ../IMET1v2_bt2/TSSUP.primary.bed")
    exit(0)

CX_report = sys.argv[1]
comp_bed = sys.argv[2]

beds = dict()
with open(comp_bed) as fh:
    for rec in fh:
        cline = rec.strip().split("\t")
        if cline[0] not in beds.keys():
            beds[cline[0]] = []
        beds[cline[0]].append(cline)

recs = []
with open(CX_report) as fh:
    for rec in fh:
        cline = rec.strip().split("\t")
        if (int(cline[3])+int(cline[4])>0):
            recs.append(cline)

def calc(rec):
    MLmc_ss = -1
    MLmc_nss = -1
    if rec[0] not in beds:
        return (MLmc_ss, MLmc_nss)
    target = beds[rec[0]]
    for fragment in target:
        if rec[2]==fragment[5] and int(rec[1])>=int(fragment[1]) and int(rec[1])<int(fragment[2]): # find match, calc it, only first occurence
            MLmc_ss = int(rec[3])*1.0/(int(rec[3])+int(rec[4]))
        if int(rec[1])>=int(fragment[1]) and int(rec[1])<int(fragment[2]):
            MLmc_nss = int(rec[3])*1.0/(int(rec[3])+int(rec[4]))
 
    return (MLmc_ss, MLmc_nss)

pool = Pool(48)
results = pool.map(calc, recs)

MLmc_gc_ss = 0
Num_CGs_gc_ss = 0
MLmc_gc_nss = 0
Num_CGs_gc_nss = 0
for mlmc_ss,mlmc_nss in results:
    if mlmc_ss>=0:
        MLmc_gc_ss += mlmc_ss
        Num_CGs_gc_ss += 1
    if mlmc_nss>=0:
        MLmc_gc_nss += mlmc_nss
        Num_CGs_gc_nss += 1

print("%.5f\t%.5f\n"%(MLmc_gc_ss*100.0/Num_CGs_gc_ss, MLmc_gc_nss*100.0/Num_CGs_gc_nss))

