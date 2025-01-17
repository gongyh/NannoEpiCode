#!/usr/bin/env python

## calc methylation level of a class of genomic component

from __future__ import print_function
import sys
from multiprocessing import Pool
import pandas as pd

if len(sys.argv) != 5:
    print("Usage: python calcMLgc.py ../IMET1v2_bt2/TSSUP.primary.bed IMET1v2.24_3_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt ...")
    exit(0)

comp_bed = sys.argv[1]
CX_report1 = sys.argv[2]
CX_report2 = sys.argv[3]
CX_report3 = sys.argv[4]

beds = []
with open(comp_bed) as fh:
    for rec in fh:
        cline = rec.strip().split("\t")
        beds.append(cline)

#chr6    0       +       0       0       CHH     CCC
CX1 = pd.read_csv(CX_report1, sep='\t', header=None, names=["Chr","Pos","Strand","mCs","Cs","Context","Seq"])
CX2 = pd.read_csv(CX_report2, sep='\t', header=None, names=["Chr","Pos","Strand","mCs","Cs","Context","Seq"])
CX3 = pd.read_csv(CX_report3, sep='\t', header=None, names=["Chr","Pos","Strand","mCs","Cs","Context","Seq"])

#CX1.describe()
#CX2.describe()
#CX3.describe()

def calc(gene):
    chrom = gene[0]
    start = int(gene[1])
    end = int(gene[2])
    name = gene[3]
    strand = gene[5]
    if strand == "+":
        end = start
        start -= 1000
    else:
        start = end
        end += 1000
    CX1_gene = CX1[(CX1["Chr"]==chrom)&(CX1["Pos"]>=start)&(CX1["Pos"]<end)]
    CX2_gene = CX2[(CX2["Chr"]==chrom)&(CX2["Pos"]>=start)&(CX2["Pos"]<end)]
    CX3_gene = CX3[(CX3["Chr"]==chrom)&(CX3["Pos"]>=start)&(CX3["Pos"]<end)]
    assert len(CX1_gene) == len(CX2_gene)
    assert len(CX3_gene) == len(CX2_gene)
    #print(CX1_gene)
    CX=pd.concat([CX1_gene,CX2_gene,CX3_gene])
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

