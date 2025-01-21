#!/usr/bin/env python3.6

## split genes bed according to DEG: up, down and no

import sys

all_genes = sys.argv[1] # genes2.bed
DE = sys.argv[2] # H0_H24_DEG.txt
prefix = sys.argv[3] # DE

all_bed = {}
with open(all_genes) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        gid = cline[3].split(".")[0]
        all_bed[gid] = line

with open(DE) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        gid = cline[0]
        if not gid.startswith("NO"):
            continue
        de = cline[1]
        fn = prefix + "_" + de + ".bed"
        with open(fn,"a") as fh2:
            fh2.write(all_bed[gid])

