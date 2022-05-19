#!/home/caoy7/anaconda2/envs/maestlin/bin/python
#--coding:utf-8 --
"""
getLocalIDS.py
Caculate the interaction density score with small window size as 1kb.


"""
__date__ = "2019-09-12"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os
import argparse
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
    Create the command line interface for the script of calcLocalIDS.py
    """
    description = """
        Caculate the interaction density score for a specific regions. 
        The output .bdg is the bedGraph result for the regions with insulation score.
        Example:
        getLocalIDS.py -f GM12878_Trac/chr21-chr21.ixy -o GM12878_Trac_chr21
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
        default=1000,
        type=int,
        help="Bin size (bp) to generate the contact matrix, default is 1000 bp."
    )
    parser.add_argument(
        "-ext",
        dest="ext",
        required=False,
        default=10,
        type=int,
        help=
        "The extension fold of the target region to show the interactions around. Default is 10."
    )
    parser.add_argument(
        "-cut",
        dest="cut",
        required=False,
        default=0,
        type=int,
        help=
        "Filtering PETs with distance < cut. Default is 0 without filtering.")

    op = parser.parse_args()
    return op


def calcIDS(f, fout, start=-1, end=-1, bs=1000, ext=10, cut=0):
    """
    Caculate the interaction density score for a region. 
    """
    print("loading %s" % f)
    key, mat = parseIxy(f, cut=cut)
    xy = XY(mat[:, 0], mat[:, 1])
    if key[0] != key[1]:
        print(
            "IS can be only caculated for intra-chromosomal interactions. Return."
        )
        return
    if start == -1:
        start = np.min(xy.xs) + bs * ext / 2
    if end == -1:
        end = np.max(xy.ys) - bs * ext / 2
    bins = int((end - start) / bs)

    print("caculating from %s to %s of %s bins" % (start, end, bins))
    with open(fout + "_ids.bdg", "w") as fo:
        for i in tqdm(range(bins)):
            x = start + i * bs
            y = start + (i + 1) * bs
            r = 0
            for j in range(int(-ext / 2), int(ext / 2 + 1)):
                if j == 0:
                    continue
                s = x + j * bs
                e = y + j * bs
                ra, rb, rab = xy.queryLoop(x, y, s, e)
                #print(x,y,s,e,rab)
                r += len(rab)
            r = r / xy.number * 10**6
            line = [key[0], x, y, r]
            fo.write("\t".join(list(map(str, line))) + "\n")


def main():
    op = help()
    calcIDS(op.fixy,
            op.output,
            start=op.start,
            end=op.end,
            bs=op.binSize,
            ext=op.ext,
            cut=op.cut)


if __name__ == "__main__":
    main()
