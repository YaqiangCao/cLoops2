#!/home/caoy7/anaconda3/envs/astroBoy/bin/python
#--coding:utf-8 --
"""
getDI.py
Calculating the Directionality Index according to the paper of Topological Domains in Mammalian Genomes Identified by Analysis of Chromatin Interactions
DI =  (B-A) / (|B-A|) * ( (A-E)**2 + (B-E)**2 ) / E
E = (A+B)/2
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
        Caculate the Directionality Index for a specific region. 
        The output .bdg is the bedGraph result for the regions with Directionality Index.
        DI is defined as accroding to the formula in following paper:
        Topological Domains in Mammalian Genomes Identified by Analysis of Chromatin Interactions
        Example:
        getDI.py -f GM12878_Trac/chr21-chr21.ixy -o GM12878_Trac_chr21
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


def calcDI(f, fout, start=-1, end=-1, bs=10000, step=100000):
    """
    Calculation of insulation score, output as .bedGraph file.
    """
    print("loading %s" % f)
    key, mat = parseIxy(f, cut=0)
    xy = XY(mat[:, 0], mat[:, 1])
    if key[0] != key[1]:
        print(
            "DI can be only caculated for intra-chromosomal interactions. Return."
        )
        return
    if start == -1:
        start = np.min(xy.xs) + step
    if end == -1:
        end = np.max(xy.ys) - step
    bins = int((end - start) / bs)
    print("caculating from %s to %s of %s bins" % (start, end, bins))
    ds = []
    ss = []
    with open(fout + ".bdg", "w") as fo:
        for i in tqdm(range(bins)):
            x = start + i * bs
            a = len(xy.queryPeakBoth(x - step, x))
            b = len(xy.queryPeakBoth(x, x + step))
            e = (a + b) / 2
            if e == 0:
                continue
            if a == b:
                di = 0
            else:
                di = (b - a) / np.abs(b - a) * ((a - e)**2 + (b - e)**2) / e
            line = [key[0], x, x + bs, di]
            fo.write("\t".join(list(map(str, line))) + "\n")


def main():
    op = help()
    calcDI(op.fixy,
           op.output,
           start=op.start,
           end=op.end,
           bs=op.binSize,
           step=op.step)


if __name__ == "__main__":
    main()
