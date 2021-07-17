#!python
#--coding:utf-8 --
"""
getSigEnrich.py
Get the enrichment of interaction signals, just like the fingerprint plot for the ChIP-seq.
"""
__date__ = "2019-08-26"
__modified__ = ""
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
    Create the command line interface for the script of estSigEnrich.py.
    """
    description = """
        Get the observed and expected enrichment trend plot based on contact matrix.
        Example:
        getSigEnrich.py -d GM12878_Trac -o GM12878_Trac -cut 0 -p 10 
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
        "Bin size (bp) to generate the contact matrix for estimation, default is 5000 bp."
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


def getObsPETs(mat, binSize=5000):
    """
    Get the number of PETs in bins.
    @param mat: [[x,y]]
    @param binSize:int, contact matrix bin size
    """
    minC = np.min(mat)
    a = (mat[:, 0] - minC) / binSize
    b = (mat[:, 1] - minC) / binSize
    a = a.astype(int)
    b = b.astype(int)
    ss = {}
    for i in range(len(a)):
        x = a[i]
        y = b[i]
        if x not in ss:
            ss[x] = {}
        if y not in ss[x]:
            ss[x][y] = 0
        ss[x][y] += 1
    sso = []
    for x in ss.keys():
        for y in ss[x].keys():
            sso.append(ss[x][y])
    return sso


def getSortBins(ds, bins=100):
    """
    Furthur bin the signal in contact matrix into bins, only care of the cumutative trend.
    """
    #default is ascending sort
    ds = np.sort(ds)
    #bin the contacts into 100 bins for comparing signal enrichment between samples
    nn = []
    step = int(len(ds) / bins)
    for i in range(0, len(ds), step):
        if i + step > len(ds):
            break
        nn.append(ds[i:i + step].sum())
    nn = np.array(nn)
    nn = np.cumsum(nn) / float(nn.sum())
    return nn


def preObs(f, cut=0, binSize=5000):
    chrom, mat = parseIxy(f, cut=cut)
    return getObsPETs(mat, binSize=binSize)


def preExp(f, cut=0, binSize=5000):
    chrom, mat = parseIxy(f, cut=cut)
    #shuffle data
    a = mat[:, 0]
    b = mat[:, 1]
    np.random.shuffle(a)
    np.random.shuffle(b)
    mat[:, 0] = a
    mat[:, 1] = b
    return getObsPETs(mat, binSize=binSize)


def plotObsExpSigEnrichment(f):
    """
    Plot the signal enrichment.
    """
    mat = pd.read_csv(f, sep="\t", index_col=0)
    fig, ax = pylab.subplots()
    ax.plot(mat.index,
            mat["observed"] * 100,
            color=colors[0],
            label="observed")
    ax.plot(mat.index,
            mat["expected"] * 100,
            color=colors[1],
            label="expected")
    ax.set_xlabel("Percentage of Bins")
    ax.set_ylabel("Percetange of PETs")
    ax.legend(loc="upper left")
    #ax2 = ax.twinx()
    #ax2.plot( mat.index, mat["Obs/Exp"],color=colors[2],label="Obs/Exp")
    #ax2.set_ylabel("Obs/Exp")
    #for t in ax2.get_yticklabels():
    #    t.set_color(colors[2])
    pylab.savefig("%s.pdf" % (f.replace(".txt", "")))


def main():
    op = help()
    if op.chroms == "":
        chroms = []
    else:
        chroms = set(op.chroms.split(","))
    fs = glob("%s/*.ixy" % op.dir)
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

    print("%s \t Getting the observed enrichment trend for %s." %
          (datetime.now(), op.dir))
    ds = Parallel(n_jobs=op.cpu,backend="multiprocessing")(
        delayed(preObs)(f, cut=op.cut, binSize=op.binSize) for f in tqdm(fs))
    ds = np.concatenate(ds)
    nn = getSortBins(ds)
    del ds

    print("%s \t Getting the expected enrichment trend for %s." %
          (datetime.now(), op.dir))
    for i in tqdm(range(op.repeats)):
        ds = Parallel(n_jobs=op.cpu,backend="multiprocessing")(
            delayed(preExp)(f, cut=op.cut, binSize=op.binSize) for f in fs)
        ds = np.concatenate(ds)
        nni = getSortBins(ds)
        del ds
        if i == 0:
            nnp = nni
        else:
            nnp += nni
    nnp = nnp / op.repeats

    ds = pd.DataFrame({"observed": nn, "expected": nni})
    #ds["Obs/Exp"] = ds["observed"]/ds["expected"]
    ds.to_csv("%s_sigEnrich.txt" % op.output, sep="\t", index_label="bins")

    if op.plot:
        plotObsExpSigEnrichment(op.output + "_sigEnrich.txt")


if __name__ == "__main__":
    main()
