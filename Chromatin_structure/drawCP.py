#! /usr/bin/env python

## draw distance based contact probability curves
## originally forked form sparseToDense.py

import sys
import numpy as np
import os
from iced import io
from iced.utils import extract_sub_contact_map
from scipy import sparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from iced.utils import get_inter_mask

def load_bed(filename):
    '''
    Load a BED file using numpy

    Parameters
    ----------
    filename : str,
        path to the file to load. The file should be a bed file

    Returns
    ------
    data : resolution

    '''
    data = np.genfromtxt(filename, dtype='str')
    return int(data[0,2])-int(data[0,1])

def load_lengths_perchr(filename, add_name=True):
    """
    Fast loading of the bed files

    Parameters
    ----------
    filename : str,
        path to the file to load. The file should be a bed file

    Returns
    -------
    lengths : the lengths of each chromosomes
    """
    data = np.loadtxt(filename, dtype="str")
    u, idx = np.unique(data[:, 0], return_index=True)
    lengths = [(data[:, 0] == i).sum() for i in u[np.argsort(idx)]]
    if add_name:
        return (np.array(lengths), u[np.argsort(idx)])
    else:
        return np.array(lengths)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("-b", "--bins", help="BED file with bins coordinates. If provided the chromosome lengths are used to define the output matrix size")
    parser.add_argument("-o", "--output", help="Output filename")

    args = parser.parse_args()

    if args.bins is None:
        print "--bins parameter is required"
        sys.exit(1)
    
    chr_lengths = load_lengths_perchr(args.bins)
    lengths = chr_lengths[0]
    chrnames = chr_lengths[1]

    ## Load counts in sparse format
    counts = io.load_counts(args.filename, lengths=lengths)
    
    ## Genome-wide dense matrix
    counts = counts.toarray()
    counts = counts + counts.T
    counts[np.diag_indices_from(counts)] /= 2
    counts = np.round(counts, 3)

    ## Output name for saving densed matrix
    output_name = os.path.basename(args.filename)
    output_name = output_name.replace(".matrix", "_dense.matrix")

    ## Save
    #np.savetxt(output_name, counts, '%s', delimiter="\t")

    resolution = load_bed(args.bins)

    ## extract every sub contact and draw
    pp = PdfPages(args.output)
    for i in range(len(chrnames)):
        chrname = chrnames[i]
        scount, slength = extract_sub_contact_map(counts, lengths, [i])
        chrlen = slength[0]
        ssum = np.sum(scount)
        sprob = np.zeros(chrlen)
        #print "\n"+chrname
        for j in range(chrlen):
            tot = 0
            num = 0
            for k in range(chrlen-j):
                num += 1
                if j == 0:
                    tot += scount[k][k+j]
                else:
                    tot += scount[k][k+j]*2
            sprob[j] = tot/num/ssum
        #print sprob
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        x = np.arange(0,chrlen*resolution,resolution)
        ax.plot(x, sprob)
        plt.xscale('log',nonposx='mask')
        plt.yscale('log',nonposy='mask')
        ax.grid(True)
        ax.set_xlabel("distance (bp)")
        ax.set_ylabel("contact probability")
        ax.set_title(chrname, fontsize='large')
        pp.savefig(fig)
    pp.close()





