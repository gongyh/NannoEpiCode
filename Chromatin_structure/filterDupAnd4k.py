#!/bin/env python

## filter duplicate reads and reads with insert size less than 4kbp

import os,sys
import pysam

inf = sys.argv[1]
outf = sys.argv[2]

samfile = pysam.AlignmentFile(inf, "rb")
outfile = pysam.AlignmentFile(outf, "wb", template=samfile)

for read in samfile.fetch(until_eof=True):
    if read.is_duplicate:
        pass
    elif (read.reference_name == read.next_reference_name) and (abs(read.reference_start - read.next_reference_start) < 4000):
        pass
    else:
        outfile.write(read)

samfile.close()

