#!/usr/bin/env python
#--coding:utf-8 --
"""
getSigDist.py
check interaction signal distribution.
"""
__date__ = "2020-01-08"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os
import argparse
from glob import glob
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#cLoops
from cLoops2.io import parseIxy
from cLoops2.settings import *


def help():
    """
    Create the command line interface.
    """
    description = """
        Get the observed/expected interaction signal distribution in contact matrix. 
        Example:
        getSigDist.py -d GM12878_Trac -o GM12878_Trac -cut 0 -p 10 
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
        default=1000,
        type=int,
        help=
        "Bin size (bp) to generate the contact matrix for estimation, default is 1000 bp."
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
        default=0,
        type=int,
        help=
        "The reapet times to shuffle PETs to get the mean expected background,default is 0, set larger than 1 to get the expected result."
    )
    parser.add_argument('-plot',
                        dest="plot",
                        required=False,
                        action="store_true",
                        help="Set to plot the result.")
    parser.add_argument(
        '-log',
        dest="log",
        required=False,
        action="store_true",
        help=
        "Whether log transform the PETs in bins for plotting, set to transform."
    )
    op = parser.parse_args()
    return op


def getObsPETs(mat, binSize=1000):
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


def preObs(f, cut=0, binSize=1000):
    chrom, mat = parseIxy(f, cut=cut)
    return getObsPETs(mat, binSize=binSize)


def preExp(f, cut=0, binSize=1000):
    chrom, mat = parseIxy(f, cut=cut)
    #shuffle data
    a = mat[:, 0]
    b = mat[:, 1]
    np.random.shuffle(a)
    np.random.shuffle(b)
    mat[:, 0] = a
    mat[:, 1] = b
    return getObsPETs(mat, binSize=binSize)


def plotObsExpSigDist(so, fout, binSize=1000, cut=0, se=None, log=False):
    """
    Plot the signal enrichment.
    """
    if log:
        so = np.log10(so)
        if se is not None:
            se = np.log10(se)
    fig, ax = pylab.subplots()
    sns.kdeplot(so, color=colors[0], label="observed", ax=ax)
    if se is not None:
        sns.kdeplot(se, color=colors[1], label="expected", ax=ax)
    if log:
        xlabel = "log10(PETs) in bins"
    else:
        xlabel = "PETs in bins"
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    ax.legend(loc="upper left")
    ax.set_title("%s resolution contact matrix with PETs distance > %s" %
                 (binSize, cut))
    pylab.savefig("%s.pdf" % (fout + "_sigDist"))


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
        if len(chroms) == 0:
            nfs.append(f)
        if len(chroms) > 0 and n[0] in chroms and n[1] in chroms:
            nfs.append(f)
    fs = nfs

    print("%s \t Getting the observed signal distribution for %s." %
          (datetime.now(), op.dir))
    ds = Parallel(n_jobs=op.cpu,backend="multiprocessing")(
        delayed(preObs)(f, cut=op.cut, binSize=op.binSize) for f in tqdm(fs))
    ds = np.concatenate(ds)
    with open(op.output + "_obs_sigDis.txt", "w") as fo:
        fo.write(",".join(list(map(str, ds))) + "\n")
    if op.repeats > 0:
        print("%s \t Getting the expected signal distribution for %s." %
              (datetime.now(), op.dir))
        for i in tqdm(range(op.repeats)):
            eds = Parallel(n_jobs=op.cpu,backend="multiprocessing")(
                delayed(preExp)(f, cut=op.cut, binSize=op.binSize) for f in fs)
            eds = list(np.concatenate(eds))
            if i == 0:
                nnp = eds
            else:
                nnp.extend(eds)
        with open(op.output + "_exp_sigDis.txt", "w") as fo:
            fo.write(",".join(list(map(str, nnp))) + "\n")
    else:
        nnp = None

    if op.plot:
        plotObsExpSigDist(ds,
                          op.output,
                          binSize=op.binSize,
                          cut=op.cut,
                          se=nnp,
                          log=op.log)


if __name__ == "__main__":
    main()
