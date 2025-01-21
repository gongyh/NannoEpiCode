#!/usr/bin/env python3.6

import sys
import numpy as np

H0_5_gene_state = sys.argv[1]
H24_5_gene_state = sys.argv[2]
resultf = sys.argv[3]

H0_5_gene_state_dict = {}
with open(H0_5_gene_state) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        H0_5_gene_state_dict[cline[1]] = int(cline[0])

H24_5_gene_state_dict = {}
with open(H24_5_gene_state) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        H24_5_gene_state_dict[cline[1]] = int(cline[0])

cs_cs = np.zeros((5,5),dtype=np.int)

for k,v in H0_5_gene_state_dict.items():
    v2 = H24_5_gene_state_dict[k]
    cs_cs[v-1][v2-1] += 1
    print("%s\t%d_%d"%(k,v,v2))

#print(cs_cs)
np.savetxt(resultf, cs_cs, fmt="%d", delimiter="\t")

