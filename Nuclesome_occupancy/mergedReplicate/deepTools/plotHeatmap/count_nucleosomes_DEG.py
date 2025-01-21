#!/usr/bin/env python3

summits_HC = "H0_summit_sorted.bed"
summits_LC = "H24_summit_sorted.bed"

g1 = "DEG_Up.bed"
g2 = "DEG_Down.bed"
g3 = "DEG_Not.bed"

interval = 2000
step = 10

summit_dict1 = dict()
with open(summits_HC) as fh:
    for line in fh:
        cl = line.split("\t")
        if cl[0] not in summit_dict1.keys():
            summit_dict1[cl[0]] = [int(cl[1])]
        else:
            summit_dict1[cl[0]].append(int(cl[1]))

summit_dict2 = dict()
with open(summits_LC) as fh:
    for line in fh:
        cl = line.split("\t")
        if cl[0] not in summit_dict2.keys():
            summit_dict2[cl[0]] = [int(cl[1])]
        else:
            summit_dict2[cl[0]].append(int(cl[1]))

def count_summits(start, end, summits, interval=1000, step=10):
    if start > end:
        step = -step
    counts = []
    for i in range(start,end,step):
        frag_start = i-500
        frag_end = i+500
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
        all_summits1 = summit_dict1[chrom]
        all_summits2 = summit_dict2[chrom]
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
        counts1 = count_summits(start, end, all_summits1)
        counts2 = count_summits(start, end, all_summits2)
        counts = []
        for i in range(len(counts1)):
            counts.append(counts2[i]-counts1[i])
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
        all_summits1 = summit_dict1[chrom]
        all_summits2 = summit_dict2[chrom]
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
        counts1 = count_summits(start, end, all_summits1)
        counts2 = count_summits(start, end, all_summits2)
        counts = []
        for i in range(len(counts1)):
            counts.append(counts2[i]-counts1[i])
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
        all_summits1 = summit_dict1[chrom]
        all_summits2 = summit_dict2[chrom]
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
        counts1 = count_summits(start, end, all_summits1)
        counts2 = count_summits(start, end, all_summits2)
        counts = []
        for i in range(len(counts1)):
            counts.append(counts2[i]-counts1[i])
        g3_counts = [i+j for i,j in zip(g3_counts,counts)]
        g3_num += 1


print("Up\tDown\tNot")
for i in range(len(g1_counts)):
    print(g1_counts[i]/g1_num,g2_counts[i]/g2_num,g3_counts[i]/g3_num,sep="\t")

