#!/usr/bin/env python3

summits_file = "H0_summit_sorted.bed"

g1 = "HC2_G1.bed"
#g1 = "test.bed"
g2 = "HC2_G2.bed"
g3 = "HC2_G3.bed"
g4 = "HC2_G4.bed"

interval = 2000
step = 10

summit_dict = dict()
with open(summits_file) as fh:
    for line in fh:
        cl = line.split("\t")
        if cl[0] not in summit_dict.keys():
            summit_dict[cl[0]] = [int(cl[1])]
        else:
            summit_dict[cl[0]].append(int(cl[1]))

def count_summits(start, end, summits, interval=1000, step=10):
    if start > end:
        step = -step
    counts = []
    for i in range(start,end,step):
        frag_start = i-1000
        frag_end = i+1000
        count = 0
        for st in summits:
            if st>frag_start and st<frag_end:
                count += 1
        counts.append(count)
    return counts

g1_counts = [0]*400
g1_num = 0
with open(g1) as fh:
    for line in fh:
        cl = line.split("\t")
        chrom = cl[0]
        if not chrom.startswith('chr'):
            continue
        all_summits = summit_dict[chrom]
        start = -1
        end = -1
        if cl[5] == "-":
            TSS = int(cl[2])
            start = TSS+2000
            end = TSS-2000
        else:
            TSS = int(cl[1])
            start = TSS-2000
            end = TSS+2000
        counts = count_summits(start, end, all_summits)
        #for i in range(len(counts)):
        #    print(counts[i])
        g1_counts = [i+j for i,j in zip(g1_counts,counts)]
        g1_num += 1


g2_counts = [0]*400
g2_num = 0
with open(g2) as fh:
    for line in fh:
        cl = line.split("\t")
        chrom = cl[0]
        if not chrom.startswith('chr'):
            continue
        all_summits = summit_dict[chrom]
        start = -1
        end = -1
        if cl[5] == "-":
            TSS = int(cl[2])
            start = TSS+2000
            end = TSS-2000
        else:
            TSS = int(cl[1])
            start = TSS-2000
            end = TSS+2000
        counts = count_summits(start, end, all_summits)
        #for i in range(len(counts)):
        #    print(counts[i])
        g2_counts = [i+j for i,j in zip(g2_counts,counts)]
        g2_num += 1


g3_counts = [0]*400
g3_num = 0
with open(g3) as fh:
    for line in fh:
        cl = line.split("\t")
        chrom = cl[0]
        if not chrom.startswith('chr'):
            continue
        all_summits = summit_dict[chrom]
        start = -1
        end = -1
        if cl[5] == "-":
            TSS = int(cl[2])
            start = TSS+2000
            end = TSS-2000
        else:
            TSS = int(cl[1])
            start = TSS-2000
            end = TSS+2000
        counts = count_summits(start, end, all_summits)
        #for i in range(len(counts)):
        #    print(counts[i])
        g3_counts = [i+j for i,j in zip(g3_counts,counts)]
        g3_num += 1


g4_counts = [0]*400
g4_num = 0
with open(g4) as fh:
    for line in fh:
        cl = line.split("\t")
        chrom = cl[0]
        if not chrom.startswith('chr'):
            continue
        all_summits = summit_dict[chrom]
        start = -1
        end = -1
        if cl[5] == "-":
            TSS = int(cl[2])
            start = TSS+2000
            end = TSS-2000
        else:
            TSS = int(cl[1])
            start = TSS-2000
            end = TSS+2000
        counts = count_summits(start, end, all_summits)
        #for i in range(len(counts)):
        #    print(counts[i])
        g4_counts = [i+j for i,j in zip(g4_counts,counts)]
        g4_num += 1


print("G1\tG2\tG3\tG4")
for i in range(len(g1_counts)):
    print(g1_counts[i]/g1_num,g2_counts[i]/g2_num,g3_counts[i]/g3_num,g4_counts[i]/g4_num,sep="\t")

