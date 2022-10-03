#!/home/caoy7/anaconda2/envs/pipi/bin/python
#--coding:utf-8 --
"""
"""
__date__ = "2021-02-10"
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
from scipy import sparse
from joblib import Parallel, delayed
from matplotlib.ticker import AutoLocator

#cLoops2
from cLoops2.io import parseIxy
from cLoops2.cmat import getObsMat, getExpMat
from cLoops2.settings import *


def help():
    """
    Create the command line interface.
    """
    description = """
        Plot the difference of matrix heatmaps for the 3D genome data for two sets.
        Example:
        plotDiffHeatmap.py -fa Trac1/chr21-chr21.ixy -fb Trac2/chr21-chr21.ixy -o Trac1vs2 -cut 1000 -bs 2000 -log
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-fa",
        dest="faixy",
        required=True,
        type=str,
        help="Input .ixy file generated by cLoops2 for first file.")
    parser.add_argument(
        "-fb",
        dest="fbixy",
        required=True,
        type=str,
        help="Input .ixy file generated by cLoops2 for second file.")
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    parser.add_argument(
        "-bs",
        dest="binSize",
        required=False,
        default=5000,
        type=int,
        help=
        "Bin size/matrix resolution (bp) to generate the contact matrix for estimation, default is 5000 bp."
    )
    parser.add_argument(
        "-start",
        dest="start",
        required=False,
        type=int,
        default=0,
        help="Start genomic coordinate for the target region,default is 0.")
    parser.add_argument(
        "-end",
        dest="end",
        required=False,
        type=int,
        default=-1,
        help=
        "End genomic coordinate for the target region,default will be inferred from the data."
    )
    parser.add_argument("-log",
                        dest="log",
                        required=False,
                        action="store_true",
                        default=False,
                        help="Whether to log transform the matrix.")
    parser.add_argument(
        "-cut",
        dest="cut",
        type=int,
        default=0,
        help="PETs with distance > cut will be kept, default is 0.")
    parser.add_argument(
        "-mcut",
        dest="mcut",
        type=int,
        default=-1,
        help="PETs with distance < mcut will be kept, default is -1 no limit.")
    parser.add_argument(
        "-na",
        dest="na",
        type=str,
        default="",
        help="Sample A name, default is inferred from data directory name.")
    parser.add_argument(
        "-nb",
        dest="nb",
        type=str,
        default="",
        help="Sample B name, default is inferred from data directory name.")
    parser.add_argument(
        "-vmin",
        dest="vmin",
        type=float,
        default=None,
        help="The minimum value shown in the heatmap and colorbar.")
    parser.add_argument(
        "-vmax",
        dest="vmax",
        type=float,
        default=None,
        help="The maxmum value shown in the heatmap and colorbar.")

    op = parser.parse_args()
    return op


def getData(f, cut=0, mcut=-1,start=0, end=-1):
    """
    """
    chrom, xy = parseIxy(f, cut=cut, mcut=mcut)
    if start == 0:
        start = np.min(xy)
    if end == -1:
        end = np.max(xy)
    ps = np.where((xy[:, 0] >= start) & (xy[:, 1] <= end))[0]
    xy = xy[ps, ]
    n = os.path.split(f)[-2]
    p = os.path.abspath(f)
    p = os.path.dirname(p)
    metaf = os.path.join(p, "petMeta.json")
    meta = json.loads(open(metaf).read())
    tot = meta["Unique PETs"]
    return n, chrom, xy, tot


def plotDiffMatHeatmap(
        fa,
        fb,
        fo,
        start=0,
        end=-1,
        r=5000,
        cut=0,
        mcut=-1,
        na="",
        nb="",
        log=False,
        vmin=None,
        vmax=None,
):
    """
    Plot the contact matrix heatmaps for compare.
    """
    labela, chroma, xya, tota = getData(fa, cut, mcut,start,end)
    labelb, chromb, xyb, totb = getData(fb, cut, mcut,start,end)
    if chroma != chromb:
        print("ERROR! %s and %s are not the same target chromosome, return." %
              (fa, fb))
        return
    if start == 0:
        start = min(np.min(xya), np.min(xyb))
    if end == -1:
        end = max(np.max(xya), np.max(xyb))
    if na == "":
        na = labela
    if nb == "":
        nb = labelb
    mata = getObsMat(xya, start, end, r)
    matb = getObsMat(xyb, start, end, r)
    sf = tota / totb
    mata = mata / sf
    if log:
        mat = np.log2((mata + 1) / (matb + 1))
        label = "log2( %s/%s )" % (na, nb)
    else:
        mat = mata - matb
        label = "%s-%s" % (na, nb)

    hights = 4
    hr = [6, 0.1]
    fig = pylab.figure(figsize=(4, hights))
    gs = mpl.gridspec.GridSpec(len(hr),
                               1,
                               height_ratios=hr,
                               top=0.95,
                               bottom=0.05,
                               left=0.1,
                               right=0.9,
                               wspace=0.05)
    pylab.suptitle("%s-%s, %s:%s-%s" % (na, nb, chroma[0], start, end),
                   fontsize=8)
    cmap = sns.color_palette("RdBu_r", 11).as_hex()
    cmap[int(len(cmap) / 2)] = "#FFFFFF"
    cmap = ListedColormap(cmap)
    ax = fig.add_subplot(gs[-2])
    cax = fig.add_subplot(gs[-1])
    sns.set(font_scale=0.5)
    ax = sns.heatmap(mat,
                     xticklabels=False,
                     yticklabels=False,
                     linewidths=0.0,
                     square=True,
                     cmap=cmap,
                     ax=ax,
                     center=0,
                     vmin=vmin,
                     vmax=vmax,
                     cbar_ax=cax,
                     cbar_kws={
                         'label': label,
                         'orientation': 'horizontal',
                         "shrink": 0.5,
                         "fraction": 0.2,
                         "anchor": (0.0, 1.0)
                     })
    cax.tick_params(labelsize=4)
    #draw the box
    ax.axvline(x=ax.get_xlim()[0], color="k", linewidth=2)
    ax.axvline(x=ax.get_xlim()[1], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[0], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[1], color="k", linewidth=2)
    pylab.savefig(fo + "_compareMatrix.pdf")


def main():
    op = help()
    plotDiffMatHeatmap(
        op.faixy,
        op.fbixy,
        op.output,
        start=op.start,
        end=op.end,
        r=op.binSize,
        cut=op.cut,
        log=op.log,
        na=op.na,
        nb=op.nb,
        vmin=op.vmin,
        vmax=op.vmax,
    )


if __name__ == "__main__":
    main()
