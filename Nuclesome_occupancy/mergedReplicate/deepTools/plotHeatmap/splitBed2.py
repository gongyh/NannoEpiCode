#!/usr/bin/env python3.6

import sys

all_genes = sys.argv[1] # genes2.bed
TPM = sys.argv[2] # HC.TMM.TPM.matrix
prefix = sys.argv[3] # HC

all_bed = {}
with open(all_genes) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        gid = cline[3].split(".")[0]
        all_bed[gid] = line

with open(TPM) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        gid = cline[0]
        tpm = float(cline[1])
        fn = prefix + "_"
        if (tpm>50):
            fn += "G1"
        elif (tpm>10) and (tpm<=50):
            fn += "G2"
        elif (tpm>1) and (tpm<=10):
            fn += "G3"
        elif (tpm<=1):
            fn += "G4"
        fn += ".bed"
        with open(fn,"a") as fh2:
            fh2.write(all_bed[gid])

