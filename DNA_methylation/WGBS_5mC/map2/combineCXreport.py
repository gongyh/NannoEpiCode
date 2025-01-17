#!/usr/bin/env python

## combine CX_report files

from __future__ import print_function
import sys

if len(sys.argv)<=1:
    print("Usage: python combineCXreport.py IMET1v2.24_1_R1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.CX_report.txt ...")
    exit(0)

result_dict = {}
i = 0
for i in range(1,len(sys.argv)):
    cx_report = open(sys.argv[i])
    for line in cx_report:
        i += 1
        if i%100000==0:
            print("+", file=sys.stderr, end='')
        cl = line.strip().split("\t")
        index = ",".join([cl[0],cl[1],cl[2],cl[5],cl[6]])
        if index not in result_dict:
            result_dict[index] = [int(cl[3]),int(cl[4])]
        else:
            result_dict[index][0] += int(cl[3])
            result_dict[index][1] += int(cl[4])
    cx_report.close()
    print("\n-------------", file=sys.stderr)

for k,v in result_dict.items():
    kl = k.split(",")
    print(kl[0]+"\t"+kl[1]+"\t"+kl[2]+"\t%d\t%d\t"%(v[0],v[1])+kl[3]+"\t"+kl[4])

