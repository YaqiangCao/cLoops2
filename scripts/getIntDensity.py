#!/usr/bin/env python
#--coding:utf-8 --
"""
getIntDensity.py
Get the interaction density for a region.
"""

__date__ = "2019-10-08"
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
from scipy.stats import hypergeom, binom, poisson

#cLoops2
from cLoops2.ds import XY
from cLoops2.io import parseTxt2Loops, ixy2pet
from cLoops2.callCisLoops import getPerRegions, estAnchorSig
from cLoops2.settings import *


def help():
    """
    Create the command line interface for the script of getAggLoopsPlot.py.
    """
    description = """
        Get the interaction density for regions.
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-d",
                        dest="predir",
                        required=True,
                        type=str,
                        help="Directory for cLoops2 pre generated.")
    parser.add_argument(
        "-b",
        dest="fbed",
        required=True,
        type=str,
        help=
        "The .bed file which contains regions to get the interaction density.")
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    parser.add_argument(
        "-pcut",
        dest="pcut",
        type=int,
        default=0,
        help=
        "Distance cutoff for PETs to filter, default is 0. Can be set as the estimated self-ligation distance cutoff."
    )
    parser.add_argument('-p',
                        dest="cpu",
                        required=False,
                        default=1,
                        type=int,
                        help="Number of CPUs to run the job, default is 1.")
    op = parser.parse_args()
    return op


def quantifyRegions(key, rs, fixy, pcut=0, pseudo=1):
    """
    @param key: str, such as chr21-chr21
    @param loops: list of Loop object
    @param fixy: cLoops2 pre generated .ixy file
    """
    print("%s\t quantify interaction density of %s regions in %s." %
          (datetime.now(), len(rs), key))
    xy = ixy2pet(fixy, cut=pcut)
    N = xy.number
    ds = {}
    for r in tqdm(rs):
        local = xy.queryPeakBoth(int(r[1]), int(r[2]))
        a = xy.queryPeak(int(r[1]), int(r[2]))
        distal = a.difference(local)
        ds["|".join(r)] = {
            "chrom":
            r[0],
            "start":
            r[1],
            "end":
            r[2],
            "name":
            r[3],
            "allPETs":
            len(local) * 2 + len(distal),
            "localPETs":
            len(local) * 2,
            "distalPETs":
            len(distal),
            "allRPKM": (len(local) * 2 + len(distal)) /
            (int(r[2]) - int(r[1])) / N / 2 * 10**9,
            "localRPKM":
            len(local) * 2 / (int(r[2]) - int(r[1])) / N / 2 * 10**9,
            "distalRPKM":
            len(distal) * 2 / (int(r[2]) - int(r[1])) / N / 2 * 10**9,
        }
    return ds


def parseBed(f):
    regions = {}
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        key = line[0] + "-" + line[0]
        if key not in regions:
            regions[key] = []
        regions[key].append(line)
    return regions


def main():
    op = help()
    regions = parseBed(op.fbed)
    metaf = op.predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(regions.keys())))
    ds = Parallel(n_jobs=op.cpu,backend="multiprocessing")(delayed(quantifyRegions)(
        key,
        regions[key],
        meta["data"]["cis"][key]["ixy"],
        pcut=op.pcut,
    ) for key in keys)
    data = {}
    for d in ds:
        for k, v in d.items():
            data[k] = v
    data = pd.DataFrame(data).T
    data.to_csv(op.output + "_quant.txt", sep="\t", index_label="rid")


if __name__ == "__main__":
    main()
