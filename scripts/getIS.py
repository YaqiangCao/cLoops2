#!/usr/bin/env python
#--coding:utf-8 --
"""
getIS.py
cLoops2 calcIS.py caculation the insulation score according to https://www.nature.com/articles/nature20158
X(x,s) = number of contacts between any pair of elements in the interval (x − s, x + s )

#X(x,s) = - log2{(X(x,s) − X(x + s/2,s/2) − X(x − s/2, s/2))/X(x,s)/0.5}
modified as X(x,s) = (X(x,s) − X(x + s/2,s/2) − X(x − s/2, s/2))/(X(x+s/2,s/2)+X(x-s/2,s/2)
directly compare the spanning and independent reads

"""
__date__ = "2019-09-11"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os, argparse
from glob import glob
from collections import Counter
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm

#cLoops2
from cLoops2.ds import XY
from cLoops2.io import parseIxy
from cLoops2.settings import *


def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Caculate the insulation score for a specific region. 
        The output .bdg is the bedGraph result for the regions with insulation score.
        IS is defined as accroding to the formula in following paper:
        Capturing pairwise and multi-way chromosomal conformations using chromosomal walks

        Example:
        getIS.py -f GM12878_Trac/chr21-chr21.ixy -o GM12878_Trac_chr21
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-f",
        dest="fixy",
        required=True,
        type=str,
        help=
        "Input .ixy file generated by cLoops2 to caculate insulation score.")
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    parser.add_argument(
        "-start",
        dest="start",
        required=False,
        default=-1,
        type=int,
        help=
        "Start genomic coordinate for the target region, default is the minmial corrdinate found in the file."
    )
    parser.add_argument(
        "-end",
        dest="end",
        required=False,
        default=-1,
        type=int,
        help=
        "End genomic coordinate for the target region, default is the maxmial corrdinate found in the file."
    )
    parser.add_argument(
        "-bs",
        dest="binSize",
        required=False,
        default=5000,
        type=int,
        help="Bin size (bp) to generate the contact matrix, default is 5000 bp."
    )
    parser.add_argument(
        "-s",
        dest="step",
        required=False,
        default=100000,
        type=int,
        help=
        "The upstream and downstream extension to caculate insulaiton score, default is 100000 bp."
    )
    op = parser.parse_args()
    return op


def calcIS(f, fout, start=-1, end=-1, bs=10000, step=100000):
    """
    Calculation of insulation score, output as .bedGraph file.
    """
    print("loading %s" % f)
    key, mat = parseIxy(f, cut=0)
    xy = XY(mat[:, 0], mat[:, 1])
    if key[0] != key[1]:
        print(
            "IS can be only caculated for intra-chromosomal interactions. Return."
        )
        return
    if start == -1:
        start = np.min(xy.xs) + step
    if end == -1:
        end = np.max(xy.ys) - step
    bins = int((end - start) / bs)
    print("caculating from %s to %s of %s bins" % (start, end, bins))
    ss = []
    ds = []
    for i in tqdm(range(bins)):
        x = start + i * bs
        xc = len(xy.queryPeakBoth(x - step, x + step))
        if xc == 0:
            continue
        xcright = len(xy.queryPeakBoth(x, x + step))
        xcleft = len(xy.queryPeakBoth(x - step, x))
        if xcright + xcleft == 0:
            continue
        xcbridge = xc - xcright - xcleft
        s = xcbridge / (xcright + xcleft)
        line = [key[0], x, x + bs, s]
        ds.append( line )
        ss.append( s ) 
    ss = np.array(ss)
    ss = (ss - np.mean(ss))/np.std(ss)
    for i in range(len(ds)):
        ds[i][-1] = ss[i]
    with open(fout + ".bdg", "w") as fo:
        for line in ds:
            fo.write("\t".join(list(map(str, line))) + "\n")


def main():
    op = help()
    calcIS(op.fixy,
           op.output,
           start=op.start,
           end=op.end,
           bs=op.binSize,
           step=op.step)


if __name__ == "__main__":
    main()
