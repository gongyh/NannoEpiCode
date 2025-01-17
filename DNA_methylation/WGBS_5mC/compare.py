#!/usr/bin/env python3

import sys

hc = sys.argv[1]
lc = sys.argv[2]

hc_dict = {}
with open(hc) as fh:
    for line in fh:
        cl = line.strip().split("\t")
        gene = cl[0].split('.')[0]
        hc_dict[gene] = float(cl[1])

diff_dict = {}
with open(lc) as fh:
    for line in fh:
        cl = line.strip().split("\t")
        gene = cl[0].split('.')[0]
        if gene in hc_dict:
            diff_dict[gene] = float(cl[1])-hc_dict[gene]
            print("%s\t%.3f"%(gene,diff_dict[gene]))

