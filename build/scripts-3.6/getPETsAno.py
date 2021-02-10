#!python
#--coding:utf-8 --
"""
getPETsAno.py
Get the annotation of PETs for enhancer/promoter and plot the stats.
"""

__date__ = "2019-10-11"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os
import sys
import json
import argparse
from glob import glob
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#cLoops2
from cLoops2.ds import XY
from cLoops2.io import parseIxy
from cLoops2.settings import *


def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Quantify the loop density.
        Get the annotations of PETs located for enhancer-promoter, enhancer-enhancer,
        promoter-promoter inteaction ratios.
        Example:
        getPETsAno.py -d GM12878_Trac -e enhancer.bed -p promoter.bed -cut 10000 -o GM12878_Trac_PETs_ano
        """
    parser = argparse.ArgumentParser(description=description)
                                     #formatter_class=RawTextHelpFormatter)
    parser.add_argument("-d",
                        dest="predir",
                        required=True,
                        type=str,
                        help="Directory for cLoops2 pre generated.")
    parser.add_argument("-enhancer",
                        dest="fe",
                        required=True,
                        type=str,
                        help="The enhancer annotation bed file.")
    parser.add_argument("-promoter",
                        dest="fp",
                        required=True,
                        type=str,
                        help="The enhancer annotation bed file.")
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    parser.add_argument("-p",
                        dest="cpu",
                        required=False,
                        type=int,
                        default=1,
                        help="Number of CPU to run the job. Default is 1.")
    parser.add_argument(
        "-pcut",
        dest="pcut",
        type=int,
        default=0,
        help=
        "Distance cutoff for PETs to filter, default is 0. Can be set as the estimated self-ligation distance cutoff. "
    )
    op = parser.parse_args()
    return op


def buildFeature(f):
    """
    For the features in bed file, map to genomic location easy for access.
    For the genomic region, if has feature, it will be shown as True.
    """
    print("%s\tBuilding coverage features of %s." % (datetime.now(), f))
    cov = {}
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        if chrom not in cov:
            cov[chrom] = set()
        #if len( cov[chrom] ) < end:
        #    cov[chrom].extend( [False] * (end-len(cov[chrom])+1))
        #for i in range( start,end):
        #    cov[chrom][i] = True
        cov[chrom].update(range(start, end))
    return cov


def findFeature(cov, chrom, start, end):
    """
    Judge if the target region overlap with features.
    """
    if chrom not in cov:
        return False
    for t in range(start, end):
        if t in cov[chrom]:
            return True
    return False


def annotatePETs(key, fixy, ecov, pcov, cut=0, ext=50):
    """
    Annotate the PETs to enhancer, promoter or Other.
    """
    print("%s\tAnnotating PETs for %s with cut > %s." %
          (datetime.now(), key, cut))
    ep = 0  #enhancer promoter
    ee = 0  #enhancer enhancer
    pp = 0  #promoter promoter
    en = 0  #enhancer none
    pn = 0  #promoter none
    nn = 0  #none none
    key2, mat = parseIxy(fixy, cut=cut)
    if key2[0] not in ecov and key2[0] not in pcov and key2[
            1] not in ecov and key2[1] not in pcov:
        return None
    for x, y in tqdm(mat):
        fae = findFeature(ecov, key2[0], x - ext, x + ext)
        fap = findFeature(pcov, key2[0], x - ext, x + ext)
        fbe = findFeature(ecov, key2[1], y - ext, y + ext)
        fbp = findFeature(pcov, key2[1], y - ext, y + ext)
        if fae == True and fbp == True:
            ep += 1
        elif fap == True and fbe == True:
            ep += 1
        elif fae == True and fbe == True:
            ee += 1
        elif fap == True and fbp == True:
            pp += 1
        elif fae == True and fbe == False and fbp == False:
            en += 1
        elif fbe == True and fae == False and fap == False:
            en += 1
        elif fap == True and fbe == False and fbp == False:
            pn += 1
        elif fbp == True and fae == False and fap == False:
            pn += 1
        else:
            nn += 1
    ss = {
        "Enhancer-Promoter": ep,
        "Enhancer-Enhancer": ee,
        "Promoter-Promoter": pp,
        "Enhancer-Other": en,
        "Promoter-Other": pn,
        "Other-Other": nn,
    }
    return ss


def main():
    op = help()
    metaf = op.predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    ecov = buildFeature(op.fe)
    pcov = buildFeature(op.fp)
    keys = list(meta["data"]["cis"].keys())
    ds = Parallel(n_jobs=1,backend="multiprocessing")(delayed(annotatePETs)(
    #ds = Parallel(n_jobs=op.cpu,backend="multiprocessing")(delayed(annotatePETs)(
        key, meta["data"]["cis"][key]["ixy"], ecov, pcov, cut=op.pcut)
                            for key in keys)
    ss = {}
    for d in ds:
        if d is None:
            continue
        for k, v in d.items():
            if k not in ss:
                ss[k] = v
            else:
                ss[k] += v
    ss = pd.Series(ss)
    ss.to_csv("%s_PETs_anos.txt" % op.output, sep="\t")


if __name__ == "__main__":
    main()
