#!/usr/bin/env python3.6

import sys

if len(sys.argv) < 4:
    print("Usage: python3.6 compare_H3K4me2.py H24_H3K4me2_summits_annotations.txt H0_H3K4me2_summits_annotations.txt consensus/H0_H24_DEG.txt")
    exit()

treat = sys.argv[1]
control = sys.argv[2]
DE = sys.argv[3]

cutoff = 40
outlier = 300

ctl_dict = dict()
with open(control) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        ctl_dict[cline[0]] = int(cline[1])

header = True
geneDE = dict()
with open(DE) as fh:
    for line in fh:
        if header:
            header = False
            continue
        cline = line.strip().split("\t")
        geneDE[cline[0]] = cline[1]

with open(treat) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        gene = cline[0]
        if gene in ctl_dict.keys():
            diff = int(cline[1])-ctl_dict[gene]
            if diff>outlier or diff<-outlier: # remove outlier
                continue
            ann = "Close"
            if diff > cutoff:
                ann = "Far"
            elif diff < -cutoff:
                ann = "Near"
            print(gene+"\t%d\t%s\t%s"%(diff,ann,geneDE[gene]))

