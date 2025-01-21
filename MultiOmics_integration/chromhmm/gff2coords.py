#!/bin/env python3

## using gff3 annotation file to create many coords files
'''
##gff-version 3
chr1    .       gene    363     2838    .       +       .       ID=NO01G00010;Name=NO01G00010
chr1    .       mRNA    363     2838    .       +       .       ID=NO01G00010.1;Name=NO01G00010.1;Parent=NO01G00010
chr1    .       five_prime_UTR  363     1539    .       +       .       ID=NO01G00010.1.utr5p1;Parent=NO01G00010.1
chr1    .       exon    363     1583    .       +       .       ID=NO01G00010.1.exon1;Parent=NO01G00010.1
chr1    .       CDS     1540    1583    .       +       0       ID=NO01G00010.1.cds1;Parent=NO01G00010.1
chr1    .       CDS     1709    2585    .       +       1       ID=NO01G00010.1.cds2;Parent=NO01G00010.1
chr1    .       exon    1709    2838    .       +       .       ID=NO01G00010.1.exon2;Parent=NO01G00010.1
chr1    .       three_prime_UTR 2586    2838    .       +       .       ID=NO01G00010.1.utr3p1;Parent=NO01G00010.1
'''
import sys

gff3 = sys.argv[1]
lenf = sys.argv[2]

chrlen = dict()
with open(lenf) as fh:
    for line in fh:
        cline = line.strip().split("\t")
        chrlen[cline[0]] = int(cline[1])

with open(gff3) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cline = line.strip().split("\t")
        if cline[2] == "mRNA":
            start = int(cline[3]) # left
            end = int(cline[4]) # right
            if cline[6] == "+": # forward strand
                print("TSS\t%s\t%d\t%d"%(cline[0],start-1,start))
                print("TSSanchor\t%s\t%d\t%s"%(cline[0],start-1,cline[6]))
                print("TSS2kb\t%s\t%d\t%d"%(cline[0],max(0,start-1-2000),min(start+2000,chrlen[cline[0]])))
                print("TES\t%s\t%d\t%d"%(cline[0],end-1,end))
                print("TESanchor\t%s\t%d\t%s"%(cline[0],end-1,cline[6]))
            else: # reverse strand
                print("TSS\t%s\t%d\t%d"%(cline[0],end-1,end))
                print("TSSanchor\t%s\t%d\t%s"%(cline[0],end-1,cline[6]))
                print("TSS2kb\t%s\t%d\t%d"%(cline[0],max(0,end-1-2000),min(end+2000,chrlen[cline[0]])))
                print("TES\t%s\t%d\t%d"%(cline[0],start-1,start))
                print("TESanchor\t%s\t%d\t%s"%(cline[0],start-1,cline[6]))

