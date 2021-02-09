#!/usr/bin/env python
#--coding:utf-8 --
"""
cLoops2 estDis.py
Get the observed and expected background of genomic distance vs genomic interaciton frequency.
2019-08-26: updated with parallel 
2019-12-30: update output orders
2020-10-27: integrate into cLoops2 main interface
"""

__date__ = "2019-08-23"
__modified__ = "2020-10-27"
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os
import sys
import argparse
from glob import glob
from collections import Counter
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#cLoops2
from cLoops2.io import parseIxy
from cLoops2.settings import *


def help():
    """
    Create the command line interface for the script of estObsExpDisFreq.py.
    """
    description = """
        Get the observed and expected random background of the genomic distance vs interaction frequency.
        Example:
        getObsExpDisFreq.py -d GM12878_Trac -o GM12878_Trac -cut 0 -p 10
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-d",
                        dest="dir",
                        required=True,
                        type=str,
                        help="Directory for cLoops2 pre generated.")
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    parser.add_argument(
        "-c",
        dest="chroms",
        required=False,
        default="",
        type=str,
        help=
        "Whether to process limited chroms, specify it as chr1,chr2,chr3, default is not. Use this to save time for quite big data."
    )
    parser.add_argument(
        "-bs",
        dest="binSize",
        required=False,
        default=5000,
        type=int,
        help=
        "Bin size /matrix resolution (bp) to generate the contact matrix for estimation, default is 5000 bp."
    )
    parser.add_argument(
        "-cut",
        dest="cut",
        type=int,
        default=0,
        help="Distance cutoff for PETs to filter, default is 0.")
    parser.add_argument('-p',
                        dest="cpu",
                        required=False,
                        default=1,
                        type=int,
                        help="Number of CPUs to run the job, default is 1.")
    parser.add_argument(
        '-r',
        dest="repeats",
        required=False,
        default=10,
        type=int,
        help=
        "The reapet times to shuffle PETs to get the mean expected background,default is 10."
    )
    parser.add_argument('-plot',
                        dest="plot",
                        required=False,
                        action="store_true",
                        help="Set to plot the result.")
    op = parser.parse_args()
    return op


def getObsDisFreq(mat, binSize=5000):
    """
    Get the relation between genomic distance with interactions using bin size based on Numpy
    """
    minC = np.min(mat)
    a = (mat[:, 0] - minC) / binSize
    b = (mat[:, 1] - minC) / binSize
    a = a.astype(int)
    b = b.astype(int)
    c = np.abs(a - b) * binSize
    c = c.astype(int)
    sso = Counter(c)
    return sso


def preObs(f, cut=0, mcut=-1, binSize=5000):
    """
    Process observed data.
    """
    chrom, mat = parseIxy(f, cut=cut,mcut=mcut)
    sso = getObsDisFreq(mat, binSize=binSize)
    return sso


def preExp(f, cut=0,mcut=-1, binSize=5000):
    """
    Process expected data.
    """
    chrom, mat = parseIxy(f, cut=cut,mcut=mcut)
    #shuffle data
    a = mat[:, 0]
    b = mat[:, 1]
    np.random.shuffle(a)
    np.random.shuffle(b)
    #old shuffing
    mat[:, 0] = a
    mat[:, 1] = b
    #new shuffling
    #mat[:, 1] = a
    #mat[:, 0] = b
    sso = getObsDisFreq(mat, binSize=binSize)
    return sso


def updateFreq(ssa, ssb):
    """
    Update the frequency dict
    """
    for k, v in ssb.items():
        if k not in ssa:
            ssa[k] = v
        else:
            ssa[k] = ssa[k] + v
    return ssa


def combineRs(ds):
    """
    Combine a list of dict.
    """
    rs = {}
    for d in ds:
        rs = updateFreq(rs, d)
    return rs


def plotObsExpDisFreq(f):
    """
    Plot the observed interaction frequency, expected interaction frequency and the ratio.
    """
    data = pd.read_csv(f, sep="\t", index_col=0)
    fig, ax = pylab.subplots()
    x = np.log10(data.index)
    ya = np.log10(data["observed"] / data["observed"].sum())
    yb = np.log10(data["expected"] / data["expected"].sum())
    ax.scatter(x, ya, color=colors[0], s=0.5, label="observed")
    ax.scatter(x, yb, color=colors[1], s=0.5, label="expected")
    ax.legend(loc="lower left", markerscale=5)
    ax.set_xlabel("Genomic distance, log10(bp)")
    ax.set_ylabel("Normalized interaction frequency,log10")
    ax2 = ax.twinx()
    ax2.scatter(x,
                np.log2(data["Obs/Exp"]),
                color=colors[2],
                s=2,
                label="Obs/Exp")
    ax2.set_ylabel("log2(Obs/Exp)")
    for t in ax2.get_yticklabels():
        t.set_color(colors[2])
    ax2.axhline(y=0.0, linestyle="--", linewidth=2, color="gray")
    ax2.legend(loc="upper right", markerscale=5)
    pylab.tight_layout()
    pylab.savefig("%s_obsAll.pdf" % (f.replace(".txt", "")))


def plotObsExpDisFreq2(f):
    """
    Plot the observed interaction frequency, expected interaction frequency and the ratio.
    """
    data = pd.read_csv(f, sep="\t", index_col=0)
    fig, ax = pylab.subplots()
    x = np.log10(data.index)
    ya = np.log10(data["observed"] / data["observed"].sum())
    yb = np.log10(data["expected"] / data["expected"].sum())
    ax.scatter(x, np.log2(data["Obs/Exp"]), color=colors[2], s=1, label="Obs/Exp")
    ax.axhline(y=0.0, linestyle="--", linewidth=2, color="gray")
    ax.legend(loc="upper right", markerscale=5)
    ax.set_xlabel("Genomic distance, log10(bp)")
    ax.set_ylabel("Interaction frequency, log2(Obs/Exp)")
    pylab.tight_layout()
    pylab.savefig("%s_obsExp.pdf" % (f.replace(".txt", "")))


def estDis( 
        d, 
        fout, 
        bs=5000,
        cpu=1,
        cut=0, 
        mcut=-1,
        chroms=[],
        repeats=10,
        plot=False,
    ):
    """
    Estimate the distance limitation for significant interactions over background.
    """
    if chroms == "":
        chroms = []
    else:
        chroms = set(chroms.split(","))
    fs = glob("%s/*.ixy" % d)
    nfs = []
    for f in fs:
        n = f.split("/")[-1].split(".ixy")[0].split("-")
        if n[0] != n[1]:
            continue
        if len(chroms) == 0:
            nfs.append(f)
        if len(chroms) > 0 and n[0] in chroms and n[1] in chroms:
            nfs.append(f)
    fs = nfs

    print(
        "%s \t Getting the observed genomic distance vs interaction frequency for %s."
        % (datetime.now(), d))
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(
        delayed(preObs)(f, cut=cut, mcut=mcut, binSize=bs) for f in tqdm(fs))
    obsDs = combineRs(ds)

    print(
        "%s \t Getting the expected genomic distance vs interaction frequency for %s."
        % (datetime.now(), d))
    expDs = {}
    for i in tqdm(range(repeats)):
        ds = Parallel(n_jobs=cpu,backend="multiprocessing")(
            delayed(preExp)(f, cut=cut,mcut=mcut, binSize=bs) for f in fs)
        ds = combineRs(ds)
        expDs = updateFreq(expDs, ds)

    #normalize the repeat times
    for k, v in expDs.items():
        expDs[k] = v / repeats

    #prepare for output
    ds = {"observed": obsDs, "expected": expDs}
    ds = pd.DataFrame(ds)
    ds = ds.fillna(0)
    ns = []
    for t in ds.itertuples():
        if np.sum(t[1:]) < 1:
            ns.append(t[0])
    ds = ds.drop(ns)
    ds["Obs/Exp"] = ds["observed"] / ds["expected"]

    ns = ds.index
    ns = np.sort(ns)
    ds = ds.loc[ns, ]
    ds.to_csv(fout + "_obsExpDisFreq.txt",
              sep="\t",
              index_label="distance(bp)")

    if plot:
        plotObsExpDisFreq(fout + "_obsExpDisFreq.txt")
        plotObsExpDisFreq2(fout + "_obsExpDisFreq.txt")


