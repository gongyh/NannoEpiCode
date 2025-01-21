#!/bin/env python

# change HiC-pro format to TADbit .tsv format

import os,sys

if len(sys.argv)<6:
    print "Usage: python hicpro2tadbit.py hic_results/data/*/*_allValidPairs ../IMET1_MboI.bed *_r1_map_len.txt *_r2_map_len.txt *.tsv"
    exit()

allvf = sys.argv[1] # hic_results/data/*/*_allValidPairs
ref = sys.argv[2] # IMET1_MboI.bed
r1lf = sys.argv[3] # *_r1_map_len.txt
r2lf = sys.argv[4] # *_r2_map_len.txt
outf = sys.argv[5] # *.tsv

# read ref, r1lf and r2lf files for use
frag_up = dict()
frag_dn = dict()
r1l = dict()
r2l = dict()

with open(ref,"r") as fh:
    for line in fh.readlines():
        cline = line.strip().split("\t")
        frag_up[cline[3]] = cline[1]
        frag_dn[cline[3]] = cline[2]

with open(r1lf,"r") as fh:
    for line in fh.readlines():
        cline = line.strip().split("\t")
        r1l[cline[0]] = cline[1]

with open(r2lf,"r") as fh:
    for line in fh.readlines():
        cline = line.strip().split("\t")
        r2l[cline[0]] = cline[1]


outfh = open(outf,"w")

with open(allvf,"r") as fh:
    for line in fh.readlines():
        cline = line.strip().split("\t")
        # E00510:210:HCNKLCCXY:3:1204:9354:51693  chr1    88786   -       chr27   241617  -       4206    HIC_chr1_221    HIC_chr27_548   42      42
        read_id = cline[0]
        chr1 = cline[1]
        pos1 = cline[2]
        strand1 = 1 if cline[3] == "+" else 0
        length1 = r1l[read_id]
        up1 = frag_up[cline[8]]
        dn1 = frag_dn[cline[8]]
        chr2 = cline[4]
        pos2 = cline[5]
        strand2 = 1 if cline[6] == "+" else 0
        length2 = r2l[read_id]
        up2 = frag_up[cline[9]]
        dn2 = frag_dn[cline[9]]
        outc = "\t".join(str(x) for x in [read_id,chr1,pos1,strand1,length1,up1,dn1,chr2,pos2,strand2,length2,up2,dn2])
        #print outc
        outfh.write(outc+"\n")

outfh.close()


