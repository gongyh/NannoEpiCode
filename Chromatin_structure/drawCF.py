#!/bin/env python

## draw distance based contact frequency
## python drawCF.py -b 100 hic_results/data/NannoH0/NannoH0_allValidPairs

import os,sys
import numpy as np

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("-d", "--d_min", type=int, help="Minimum distance (bp) to plot")
    parser.add_argument("-D", "--d_max", type=int, help="Maximum distance (bp) to plot")
    parser.add_argument("-b", "--bins", type=int, default=100, help="Number of bins")
    parser.add_argument("-q", "--min_qv", type=int, help="Minimum qv values of alignment")
    parser.add_argument("-f", "--max_frag", type=int, help="Maximum fragment size to consider")
    parser.add_argument("-l", "--log10", type=bool, default=True, help="log10 scale")
    parser.add_argument("-o", "--output", help="Output filename")

    args = parser.parse_args()

    infh = open(args.filename)

    ffArr = np.array([,])
    frArr = np.array([,])
    rfArr = np.array([,])
    rrArr = np.array([,])

    for line in infh:
        cline = line.strip().split('\t')

        r1_chr = cline[1]
        r1_pos = int(cline[2])
        r1_str = cline[3]
        r2_chr = cline[4]
        r2_pos = int(cline[5])
        r2_str = cline[6]
        frag_size = int(cline[7])
        r1_qv = int(cline[10])
        r2_qv = int(cline[11])

        if args.max_frag is not None and frag_size < args.max_frag: # so longer fragment, filter out
            continue

        if args.min_qv is not None and r1_qv < args.min_qv or r2_qv < args.min_qv: # low mapping quality, filter out
            continue

        if r1_chr == r2_chr: # same chromosome
            if r1_str == "+" and r2_str == "+": # FF, same
                np.append(ffArr[r1_chr],abs(r2_pos-r1_pos))
            elif r1_str == "+" and r2_str == "-": # FR, inward
                np.append(frArr[r1_chr],abs(r2_pos-r1_pos))
            elif r1_str == "-" and r2_str == "+": # RF, outward
                np.append(rfArr[r1_chr],abs(r2_pos-r1_pos))
            elif r1_str == "-" and r2_str == "-": # RR, same
                np.append(rrArr[r1_chr],abs(r2_pos-r1_pos))
            else:
                print "can not be here!"

        print ffArr
    infh.close()

    #print ffArr.size, frArr.size, rfArr.size, rrArr.size

