#!/usr/bin/env python

## calc methylation level distribution around genes

from __future__ import print_function
import sys
from multiprocessing import Pool

if len(sys.argv) != 6:
    print("Usage: python calcMLdistr.py IMET1v2.H_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt LH_genes.TMM.TPM.avr.matrix ../IMET1v2_bt2/genes.primary.bed H ../IMET1v2_bt2/IMET1.bed")
    exit(0)

CX_report = sys.argv[1] # IMET1v2.H_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt
expression = sys.argv[2] # LH_genes.TMM.TPM.avr.matrix
genes_bed = sys.argv[3] # ../IMET1v2_bt2/genes.primary.bed
sample = sys.argv[4] # H, L
genome_bed = sys.argv[5]

chr_end = dict()
with open(genome_bed) as fh:
    for line in fh:
        cline = line.strip().split('\t')
        chr_end[cline[0]] = cline[2]

beds = []
with open(genes_bed) as fh:
    for rec in fh:
        cline = rec.strip().split("\t")
        if int(cline[1])>=1500 and int(cline[2])+1500<=chr_end[cline[0]]:
            beds.append(cline)

recs = dict()
with open(CX_report) as fh:
    for rec in fh:
        cline = rec.strip().split("\t")
        if cline[0] not in recs:
            recs[cline[0]] = []
        recs[cline[0]].append(cline)

exps = dict()
header = True
headers = []
col = 1 # default column
with open(expression) as fh:
    for rec in fh:
        cline = rec.rstrip().split("\t")
        if header:
            headers = cline
            index = 0
            for h in headers:
                if h==sample:
                    col = index
                    break
                index += 1
            header = False
            continue
        tpm = float(cline[col])
        gname = cline[0].split(".")[0]
        if tpm<1:
            exps[gname] = "None"
        elif tpm >=1 and tpm < 10:
            exps[gname] = "Low"
        elif tpm >=10 and tpm < 50:
            exps[gname] = "Medium"
        elif tpm >= 50:
            exps[gname] = "High"
        else:
            pass

def fragment(start,end,strand,point):
    if point < start: # before
        index = (point - (start - 1500))//100
    elif point >= start and point < end: # in gene body
        index = 15 + (point - start)*20//(end - start)
    elif point >= end: # after
        index = 15 + 20 + (point-end)//100
    else: # will not be here forever
        pass

    if strand == "+":
        return index
    elif strand == "-":
        return 49-index
    else: # . etc non-strand-specific
        return index

def calc(bed):
    #print('.',end='',file=sys.stderr)
    MLmc_ss = [0]*50
    MLmc_nss = [0]*50
    MLmc_ss_num = [0]*50
    MLmc_nss_num = [0]*50
    target = recs[bed[0]]
    gene = bed[3].split(".")[0]
    for rec in target:
        index = fragment(int(bed[1]), int(bed[2]), bed[5], int(rec[1]))
        total = int(rec[3])+int(rec[4])
        if bed[5]==rec[2] and int(rec[1])>=int(bed[1])-1500 and int(rec[1])<int(bed[2])+1500 and total>0: # find match, calc it
            MLmc_ss[index] += int(rec[3])*1.0/total
            MLmc_ss_num[index] += 1
        if int(rec[1])>=int(bed[1])-1500 and int(rec[1])<int(bed[2])+1500 and total>0:
            MLmc_nss[index] += int(rec[3])*1.0/total
            MLmc_nss_num[index] += 1

    return (exps[gene], MLmc_ss, MLmc_nss, MLmc_ss_num, MLmc_nss_num)

pool = Pool(48)
results = pool.map(calc, beds)

MLdistr_ss_None = [0]*50
MLdistr_nss_None = [0]*50
MLdistr_ss_Low = [0]*50
MLdistr_nss_Low = [0]*50
MLdistr_ss_Medium = [0]*50
MLdistr_nss_Medium = [0]*50
MLdistr_ss_High = [0]*50
MLdistr_nss_High = [0]*50

MLdistr_CG_ss_None = [0]*50
MLdistr_CG_nss_None = [0]*50
MLdistr_CG_ss_Low = [0]*50
MLdistr_CG_nss_Low = [0]*50
MLdistr_CG_ss_Medium = [0]*50
MLdistr_CG_nss_Medium = [0]*50
MLdistr_CG_ss_High = [0]*50
MLdistr_CG_nss_High = [0]*50

for exps_gene,mlmc_ss,mlmc_nss,mlmc_ss_num,mlmc_nss_num in results:
    if exps_gene=="None":
        for i in range(0,50):
            if mlmc_ss[i]>=0:
                MLdistr_ss_None[i] += mlmc_ss[i]
                MLdistr_CG_ss_None[i] += mlmc_ss_num[i]
            if mlmc_nss[i]>=0:
                MLdistr_nss_None[i] += mlmc_nss[i]
                MLdistr_CG_nss_None[i] += mlmc_nss_num[i]
    elif exps_gene=="Low":
        for i in range(0,50):
            if mlmc_ss[i]>=0:
                MLdistr_ss_Low[i] += mlmc_ss[i]
                MLdistr_CG_ss_Low[i] += mlmc_ss_num[i]
            if mlmc_nss[i]>=0:
                MLdistr_nss_Low[i] += mlmc_nss[i]
                MLdistr_CG_nss_Low[i] += mlmc_nss_num[i]
    elif exps_gene=="Medium":
        for i in range(0,50):
            if mlmc_ss[i]>=0:
                MLdistr_ss_Medium[i] += mlmc_ss[i]
                MLdistr_CG_ss_Medium[i] += mlmc_ss_num[i]
            if mlmc_nss[i]>=0:
                MLdistr_nss_Medium[i] += mlmc_nss[i]
                MLdistr_CG_nss_Medium[i] += mlmc_nss_num[i]
    elif exps_gene=="High":
        for i in range(0,50):
            if mlmc_ss[i]>=0:
                MLdistr_ss_High[i] += mlmc_ss[i]
                MLdistr_CG_ss_High[i] += mlmc_ss_num[i]
            if mlmc_nss[i]>=0:
                MLdistr_nss_High[i] += mlmc_nss[i]
                MLdistr_CG_nss_High[i] += mlmc_nss_num[i]

ss_None = [0]*50
nss_None = [0]*50
ss_Low = [0]*50
nss_Low = [0]*50
ss_Medium = [0]*50
nss_Medium = [0]*50
ss_High = [0]*50
nss_High = [0]*50

for i in range(0,50):
    ss_None[i] = MLdistr_ss_None[i]*100.0/MLdistr_CG_ss_None[i] if MLdistr_CG_ss_None[i] > 0 else 0.0
    nss_None[i] = MLdistr_nss_None[i]*100.0/MLdistr_CG_nss_None[i] if MLdistr_CG_nss_None[i] > 0 else 0.0
    ss_Low[i] = MLdistr_ss_Low[i]*100.0/MLdistr_CG_ss_Low[i] if MLdistr_CG_ss_Low[i] > 0 else 0.0
    nss_Low[i] = MLdistr_nss_Low[i]*100.0/MLdistr_CG_nss_Low[i] if MLdistr_CG_nss_Low[i] > 0 else 0.0
    ss_Medium[i] = MLdistr_ss_Medium[i]*100.0/MLdistr_CG_ss_Medium[i] if MLdistr_CG_ss_Medium[i] > 0 else 0.0
    nss_Medium[i] = MLdistr_nss_Medium[i]*100.0/MLdistr_CG_nss_Medium[i] if MLdistr_CG_nss_Medium[i] > 0 else 0.0
    ss_High[i] = MLdistr_ss_High[i]*100.0/MLdistr_CG_ss_High[i] if MLdistr_CG_ss_High[i] > 0 else 0.0
    nss_High[i] = MLdistr_nss_High[i]*100.0/MLdistr_CG_nss_High[i] if MLdistr_CG_nss_High[i] > 0 else 0.0

print("ss_None",end='')
for i in range(0,50):
    print("\t%.5f"%ss_None[i], end='')
print("\n",end='')

print("ss_Low",end='')
for i in range(0,50):
    print("\t%.5f"%ss_Low[i], end='')
print("\n",end='')

print("ss_Medium",end='')
for i in range(0,50):
    print("\t%.5f"%ss_Medium[i], end='')
print("\n",end='')

print("ss_High",end='')
for i in range(0,50):
    print("\t%.5f"%ss_High[i], end='')
print("\n",end='')

print("nss_None",end='')
for i in range(0,50):
    print("\t%.5f"%nss_None[i], end='')
print("\n",end='')

print("nss_Low",end='')
for i in range(0,50):
    print("\t%.5f"%nss_Low[i], end='')
print("\n",end='')

print("nss_Medium",end='')
for i in range(0,50):
    print("\t%.5f"%nss_Medium[i], end='')
print("\n",end='')

print("nss_High",end='')
for i in range(0,50):
    print("\t%.5f"%nss_High[i], end='')
print("\n",end='')

