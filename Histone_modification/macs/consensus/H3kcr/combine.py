#!/usr/bin/env python3.6

import sys

if len(sys.argv) < 5:
    print("Usage: python3.6 combine.py ../H0_H24_DEG.txt H3k9ac.consensus_peaks.Nearest_PromoterID.txt Up_Peaks.txt Down_Peaks.txt")
    exit(0)

deg = sys.argv[1] # ../H0_H24_DEG.txt
pid = sys.argv[2] # H3k9ac.consensus_peaks.Nearest_PromoterID.txt
up = sys.argv[3] # Up_Peaks.txt
dn = sys.argv[4] # Down_Peaks.txt

degs = dict()
header = True
with open(deg) as fh:
    for line in fh:
        if header:
            header = False
            continue
        cl = line.strip().split("\t")
        degs[cl[0]] = cl[1]

ups = []
with open(up) as fh:
    for line in fh:
        cl = line.strip()
        ups.append(cl)

dns = []
with open(dn) as fh:
    for line in fh:
        cl = line.strip()
        dns.append(cl)

print("Peak\tGene\tDEG\tDEPeak")
with open(pid) as fh:
    for line in fh:
        cl = line.strip().split("\t")
        peak = cl[0]
        gene = cl[1]
        de_gene = degs[gene]
        de_peak = "Not"
        if peak in ups:
            de_peak = "Up"
        elif peak in dns:
            de_peak = "Down"
        print("%s\t%s\t%s\t%s" % (peak, gene, de_gene, de_peak))

