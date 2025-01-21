#!/usr/bin/env python3.6

import sys

## combine cs transitions and DEG information
# python3.6 combine.py H0_H24_CS.details HC_VLC_DEG_Up.txt HC_VLC_DEG_Down.txt

cs = sys.argv[1] # H0_H24_CS.details
up = sys.argv[2] # HC_VLC_DEG_Up.txt
dn = sys.argv[3] # HC_VLC_DEG_Down.txt

up_genes = []
header = True
with open(up) as fh:
    for line in fh:
        if header:
            header = False
            continue
        cline = line.strip()
        up_genes.append(cline)


dn_genes = []
header = True
with open(dn) as fh:
    for line in fh:
        if header:
            header = False
            continue
        cline = line.strip()
        dn_genes.append(cline)

print("Gene\tCS\tTransition\tDEG")
with open(cs) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        gid = cline[0]
        change = cline[1]
        raw = cline[1].split("_")[0]
        DEG = "Not"
        if gid in up_genes:
            DEG = "Up"
        elif gid in dn_genes:
            DEG = "Down"
        else:
            DEG = "Not"
        print("%s\t%s\t%s\t%s"%(gid,raw,change,DEG))

