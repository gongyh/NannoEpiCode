#!/usr/bin/env python

## calc methylation level every 10kb

from __future__ import print_function
import sys
from multiprocessing import Pool
import pandas as pd

if len(sys.argv) != 3:
    print("Usage: python calcMLgc.py 10kb.bed IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt")
    exit(0)

comp_bed = sys.argv[1]
CX_report1 = sys.argv[2]

beds = []
with open(comp_bed) as fh:
    for rec in fh:
        cline = rec.strip().split("\t")
        beds.append(cline)

#chr6    0       +       0       0       CHH     CCC
CX1 = pd.read_csv(CX_report1, sep='\t', header=None, names=["Chr","Pos","Strand","mCs","Cs","Context","Seq"])

#CX1.describe()

def calc(gene):
    chrom = gene[0]
    start = int(gene[1])
    end = int(gene[2])
    name = gene[3]
    strand = gene[4]
    CX1_gene = CX1[(CX1["Chr"]==chrom)&(CX1["Pos"]>=start)&(CX1["Pos"]<end)]
    #print(CX1_gene)
    CX=pd.concat([CX1_gene])
    #print(name)
    #print(CX)
    num = 0
    mcpc = 0
    for idx,row in CX.iterrows():
        if row[4]!=0:
            num += 1
            mcpc += row[3]/row[4]
    if num>=10:
        return (name, mcpc*100.0/num)
    else:
        return (name, -1)

pool = Pool(48)
results = pool.map(calc, beds)

for name,value in results:
    if value != -1:
        print(name+"\t%.3f"%value)

