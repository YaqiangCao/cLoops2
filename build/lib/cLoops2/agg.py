#!/usr/bin/env python
#--coding:utf-8--
"""
Aggregated analysis related functions.
2020-02-24: basically finished. ES added for domains. Normalizaiton added for aggregated domains.
2020-02-25: Normalization for aggregated domains not work well in the heatmap, do not try again.
2020-07-07: aggregated view points analysis added.
2020-11-10: updated view points plot, no more trying of data normalization, not work well.
2020-11-11: added new aggregation of two anchors
2020-11-12: twoAnchor aggregation analysis should be modified to have directions, such as anchor as are CTCF + and anchor bs are CTCF -
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import warnings
warnings.filterwarnings("ignore")
import os
import json
from glob import glob
from copy import deepcopy
from datetime import datetime

#3rd library
import pylab
import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from scipy.stats import zscore
from matplotlib.colors import ListedColormap
from scipy.ndimage.interpolation import rotate  #rotate the domain
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#cLoops2
from cLoops2.ds import XY, Loop
from cLoops2.io import parseIxy, parseTxt2Loops
from cLoops2.cmat import get1DSigMat, getBinMean, getObsMat
from cLoops2.utils import getLogger
from cLoops2.settings import *


### aggregated peaks releated funcitons
def getAPeaks(key,
              fixy,
              peaks,
              ext=5000,
              bins=100,
              cut=0,
              mcut=-1,
              skipZeros=False):
    """
    Get the contact matrix for target regions.
    """
    key2, mat = parseIxy(fixy, cut=cut, mcut=mcut)
    if mat.shape[0] == 0:
        print(
            "No PETs found in %s maybe due to distance cutoff for PET > %s <%s."
            % (fixy, cut, mcut))
    xy = XY(mat[:, 0], mat[:, 1])
    print("Get the 1D signal array nearby target %s regions from %s." %
          (len(peaks), key))
    for p in peaks:
        center = (p[1] + p[2]) / 2
        p[1] = int(max(0, center - ext))
        p[2] = int(center + ext)
    xs = np.arange(-ext, ext, ext * 2 / bins)
    ss = get1DSigMat(xy, peaks, bins=bins, skipZeros=skipZeros)
    if ss is not None:
        ss.columns = xs
        return ss
    else:
        return None


def plotAPeaks(mat, fout, norm=False):
    """
    Plot the mean profile and heatmap
    """
    #prepare fig
    fig = pylab.figure(
        figsize=(2, 5))  #the extra 1, 0.5 for 1D signal, 0.5 for arch
    hr = [1.5, 4, 0.1]
    gs = mpl.gridspec.GridSpec(3, 1, height_ratios=hr)

    n = mat.shape[0]

    #plot 1D signal
    ax = fig.add_subplot(gs[0])
    xs = [int(float(t)) for t in list(mat.columns)]
    ss = mat.mean(axis=0)
    ax.plot(xs, ss, color=colors[1])
    ax.set_ylim([0, np.max(ss) * 1.1])
    ax.set_title("mean 1D profile", fontsize=6, y=0.98)
    ax.tick_params(axis='both', which='major', labelsize=5)
    ax.set_xlim([np.min(xs), np.max(xs)])
    ax.set_ylabel("RPM")

    #plot mean heatmap
    #process the matrix
    #sort the matrix
    mat = mat.values
    ns = np.mean(mat, axis=1)
    ns = np.argsort(ns)
    ns = ns[::-1]
    mat = mat[ns, ]
    #bin the matrix rows from thousands to 200
    if mat.shape[0] > 200:
        nmat = []
        s = int(mat.shape[0] / 200)
        for i in range(0, 200):
            nmat.append(np.mean(mat[s * i:s * (i + 1)], axis=0))
        nmat.append(np.mean(mat[s * 200:], axis=0))
        mat = np.array(nmat)
    #mat = np.log10(mat+1)
    if norm:
        mat = zscore(mat, axis=1)
        cmap = cmaps["div"]
        center = 0
        vmin = None
        ano = "zscore of RPM"
    else:
        cmap = cmaps["red"]
        vmin = 0
        center = None
        ano = "RPM"
    ax = fig.add_subplot(gs[1])
    cax = fig.add_subplot(gs[2])
    sns.heatmap(mat,
                xticklabels=False,
                yticklabels=False,
                square=False,
                ax=ax,
                cmap=cmap,
                cbar_ax=cax,
                linewidths=0,
                center=center,
                vmin=vmin,
                cbar_kws={
                    'label': ano,
                    'orientation': 'horizontal',
                    "shrink": 0.2,
                })
    ax.set_ylabel("%s regions" % n)
    ax.axvline(x=ax.get_xlim()[0], color="k", linewidth=2)
    ax.axvline(x=ax.get_xlim()[1], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[0], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[1], color="k", linewidth=2)
    ax.set_xlabel("Distance from center (bp)")
    pylab.savefig("%s_aggPeaks.pdf" % fout)


def aggPeaks(predir,
             peakf,
             output,
             logger,
             ext=5000,
             bins=100,
             cut=0,
             mcut=-1,
             cpu=1,
             skipZeros=False,
             norm=False):
    """
    Aggregation peak analysis.
    """
    #read in peaks orgainzed by chromosomes
    peaks = {}
    for line in open(peakf):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        key = line[0] + "-" + line[0]
        if key not in peaks:
            peaks[key] = []
        line[1] = int(line[1])
        line[2] = int(line[2])
        peaks[key].append(line)
    #meta data
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(peaks.keys())))
    #get the data
    ds = Parallel(n_jobs=cpu, backend="multiprocessing")(delayed(getAPeaks)(
        key,
        meta["data"]["cis"][key]["ixy"],
        peaks[key],
        ext=ext,
        bins=bins,
        cut=cut,
        mcut=mcut,
        skipZeros=skipZeros,
    ) for key in keys)
    #save the data
    ds = [d for d in ds if d is not None]
    mat = pd.concat(ds, axis=0)
    total = meta["Unique PETs"] * 2
    mat = mat / total * 10**6
    mat.to_csv(output + "_aggPeaks.txt", sep="\t")
    logger.info(
        "%s RPM normalized 1D signal around anchors saved to %s_aggPeaks.txt" %
        (mat.shape[0], output))
    #plot the data
    plotAPeaks(mat, output, norm=norm)


### aggregated loops releated funcitons
def getPerRegions(loop, win=5):
    """
    Get the nearby regions for interacting two locus, win as how many nearby, 6 is enough for interacting more than 100 regions to estimate FDR and others. The mean distance of all the permutated regions is the same to that between iva and ivb.
    @param loop: cLoops:ds:Loop 
    @param xy: .pet file parse data
    """
    # the PET id in the permutated regions
    ivas, ivbs = [], []
    ca = loop.x_start
    cb = loop.y_start
    sa = loop.x_end - loop.x_start
    sb = loop.y_end - loop.y_start
    for i in range(0 - win, win + 1):
        niva = [max([0, ca + i * sa]), max([0, ca + (i + 1) * sa])]
        nivb = [max([0, cb + i * sb]), max([0, cb + (i + 1) * sb])]
        ivas.append(niva)
        ivbs.append(nivb)
    return ivas, ivbs


def getALoops(key,
              fixy,
              loops,
              ext=10,
              cut=0,
              mcut=-1,
              oneBins=100,
              skipZeros=False,
              oneD=False):
    """
    Get the contact matrix for target loop and nearby regions for one chromosomal.
    """
    key2, mat = parseIxy(fixy, cut=cut, mcut=mcut)
    if mat.shape[0] == 0:
        print(
            "No PETs found in %s maybe due to distance cutoff for PET > %s <%s."
            % (fixy, cut, mcut))
        return None, None
    #XY object
    xy = XY(mat[:, 0], mat[:, 1])
    print(
        "Getting the contact matrix for target %s loops and nearby regions from %s."
        % (len(loops), key))

    #matrix of loops and nearby
    mats = []
    rs = []
    for loop in tqdm(loops):
        #get the loop and nearby
        mat = np.zeros([ext + 1, ext + 1])
        ivas, ivbs = getPerRegions(loop, int(ext / 2))
        nas, nbs = [], []
        for iva in ivas:
            nas.append(xy.queryPeak(iva[0], iva[1]))
        for ivb in ivbs:
            nbs.append(xy.queryPeak(ivb[0], ivb[1]))
        for i, na in enumerate(nas):
            for j, nb in enumerate(nbs):
                nrab = float(len(na.intersection(nb)))
                mat[i][j] = nrab
        if skipZeros and np.sum(mat) == 0:
            continue
        mats.append(mat)
        #prepare for 1D profile
        startx = loop.x_start - (loop.x_end - loop.x_start) * ext
        endx = loop.x_end + (loop.x_end - loop.x_start) * ext
        starty = loop.y_start - (loop.y_end - loop.y_start) * ext
        endy = loop.y_end + (loop.y_end - loop.y_start) * ext
        rs.append([loop.chromX, startx, endx])
        rs.append([loop.chromY, starty, endy])
    mats = np.array(mats)
    if oneD:
        ss = get1DSigMat(xy, rs, bins=oneBins)
    else:
        ss = None
    return mats, ss


def getBwSig(bw, loops, bins=100, ext=10, skipZeros=False):
    """
    Get the 1D sig for bigwig around loop anchors.
    """
    print("Getting signal from %s around loop anchors." % bw)
    name = bw.split("/")[-1].split(".bw")[0]
    bw = pyBigWig.open(bw)
    ys = np.zeros(bins)
    c = 0
    for loop in tqdm(loops):
        startx = loop.x_start - (loop.x_end - loop.x_start) * ext
        endx = loop.x_end + (loop.x_end - loop.x_start) * ext
        starty = loop.y_start - (loop.y_end - loop.y_start) * ext
        endy = loop.y_end + (loop.y_end - loop.y_start) * ext
        if endx - startx < bins or endy - starty < bins:
            continue
        try:
            sx = bw.values(loop.chromX, int(startx), int(endx))
            sx = np.nan_to_num(sx)
            sx = getBinMean(sx, bins)
            sy = bw.values(loop.chromY, int(starty), int(endy))
            sy = np.nan_to_num(sy)
            sy = getBinMean(sy, bins)
            if skipZeros:
                if np.sum(sx) > 0:
                    ys = ys + sx
                    c += 1
                if np.sum(sy) > 0:
                    ys = ys + sy
                    c += 1
            else:
                ys = ys + sx + sy
                c += 2
        except:
            continue
    ys = ys / c
    return name, ys


def plotALoops( 
                mat, 
                fout, 
                ss=None, 
                bwSigs=None, 
                bins=100, 
                norm=False,
                vmin=None,
                vmax=None,
            ):
    """
    Plot the mean loop heatmap.
    bwSigs: {name:signal array}
    """
    #prepare data
    es = []
    for i in range(mat.shape[0]):
        p = int(mat[i].shape[0] / 2)
        if mat[i].mean() > 0:
            nmat = deepcopy(mat[i])
            nmat[p, p] = 0
            if nmat.mean() == 0.0:
                #es.append(0.0)
                continue
            else:
                es.append(mat[i, p, p] / nmat.mean())
        else:
            es.append(0.0)
        if norm:
            if mat[i].sum() == 0.0:
                continue
            else:
                mat[i] = mat[i] / mat[i].sum()
                mat[i] = (mat[i] - mat[i].mean()) / mat[i].std()
    mat = np.mean(mat, axis=0)
    #1D signal
    if ss is not None:
        ss = np.mean(ss, axis=0)

    #prepare fig
    if bwSigs is not None and ss is not None:
        fig = pylab.figure(figsize=(2, 3 + len(bwSigs) * 0.5))
        hr = [0.5] * len(bwSigs)
        hr.extend([0.5, 2.4, 0.1])
        gs = mpl.gridspec.GridSpec(len(bwSigs) + 3, 1, height_ratios=hr)
    elif bwSigs is None and ss is not None:
        fig = pylab.figure(figsize=(2, 3))  #the extra 1, 0.5 for 1D signal,
        hr = [0.5, 2.4, 0.1]
        gs = mpl.gridspec.GridSpec(3, 1, height_ratios=hr)
    elif bwSigs is not None and ss is None:
        fig = pylab.figure(figsize=(2, 2.5 + len(bwSigs) * 0.5))
        hr = [0.5] * len(bwSigs)
        hr.extend([2.4, 0.1])
        gs = mpl.gridspec.GridSpec(len(bwSigs) + 2, 1, height_ratios=hr)
    else:
        fig = pylab.figure(figsize=(2, 2.5))  #the extra 1, 0.5 for 1D signal,
        hr = [2.4, 0.1]
        gs = mpl.gridspec.GridSpec(2, 1, height_ratios=hr)

    xs = list(range(bins))

    #plot bigWig signals
    if bwSigs is not None:
        for i, name in enumerate(list(bwSigs.keys())):
            ax = fig.add_subplot(gs[i])
            ax.plot(xs, bwSigs[name], color=colors[i], label=name)
            ax.tick_params(axis='both', which='major', labelsize=4)
            ax.set_xticklabels([])
            ax.set_xlim([np.min(xs), np.max(xs)])
            ax.set_ylim([0, np.max(bwSigs[name]) * 1.1])
            ax.legend(fontsize=5, fancybox=False, frameon=False)

    if ss is not None:
        #plot 1D signal
        if bwSigs is not None:
            ax = fig.add_subplot(gs[len(bwSigs)])
            ax.plot(xs, ss, color=colors[len(bwSigs)], label="mean 1D profile")
        else:
            ax = fig.add_subplot(gs[0])
            ax.plot(xs, ss, color=colors[0], label="mean 1D profile")
        ax.legend(fontsize=5, fancybox=False, frameon=False)
        ax.tick_params(axis='both', which='major', labelsize=4)
        ax.set_xlim([np.min(xs), np.max(xs)])
        ax.set_ylim([0, np.max(ss) * 1.1])
        ax.set_xticklabels([])
        ax.set_ylabel("RPM", fontsize=6)

    #plot mean heatmap
    ax = fig.add_subplot(gs[-2])
    cax = fig.add_subplot(gs[-1])
    if norm == False:
        cmap = sns.cubehelix_palette(light=1, as_cmap=True)
        center = None
        if vmin is None:
            vmin = 0
    else:
        cmap = sns.color_palette("RdBu_r", 11).as_hex()
        cmap[int(len(cmap) / 2)] = "#FFFFFF"
        cmap = ListedColormap(cmap)
        center = 0
    sns.heatmap(mat,
                xticklabels=False,
                yticklabels=False,
                square=True,
                ax=ax,
                cmap=cmap,
                cbar_ax=cax,
                linewidths=0.05,
                linecolor="gray",
                linestyle="--",
                vmin=vmin,
                vmax=vmax,
                center=center,
                cbar_kws={
                    'orientation': 'horizontal',
                    "shrink": 0.2,
                })
    cax.tick_params(labelsize=4)
    ax.axvline(x=ax.get_xlim()[0], color="k", linewidth=2)
    ax.axvline(x=ax.get_xlim()[1], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[0], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[1], color="k", linewidth=2)
    ax.set_title("%s loops; ES:%.3f" % (len(es), np.mean(es)),
                 fontsize=6,
                 y=0.97)
    pylab.savefig("%s_aggLoops.pdf" % fout)


def aggLoops(predir,
             loopf,
             output,
             logger,
             bws="",
             ext=10,
             cut=0,
             mcut=-1,
             lcut=0,
             cpu=1,
             skipZeros=False,
             norm=False,
             oneD=False,
             vmax=None,
             vmin=None,
             ):
    """
    Aggregated loops analysis. 
    """
    if bws != "":
        bws = bws.split(",")
        bws = [bw for bw in bws if os.path.isfile(bw)]
    else:
        bws = []

    loops = parseTxt2Loops(loopf, cut=lcut)
    if len(loops) == 0:
        logger.info("No loops found in %s with loop ditance > %s, return." %
                    (loopf, lcut))
        return

    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(loops.keys())))
    ds = Parallel(n_jobs=cpu, backend="multiprocessing")(delayed(getALoops)(
        key,
        meta["data"]["cis"][key]["ixy"],
        loops[key],
        ext=ext,
        cut=cut,
        mcut=mcut,
        skipZeros=skipZeros,
        oneD=oneD,
    ) for key in keys)
    """
    #trans
    keys = list(meta["data"]["trans"].keys())
    keys = list(set(keys).intersection(set(loops.keys())))
    ds.extend(
        Parallel(n_jobs=op.cpu,backend="multiprocessing")(
            delayed(getSubMat)(key,
                               meta["data"]["trans"][key]["ixy"],
                               loops[key],
                               ext=op.ext,
                               cut=op.pcut) for key in keys))
    """
    #matrix
    mat = np.concatenate([d[0] for d in ds if d[0] is not None], axis=0)
    np.savez(output + "_aggLoops.npz", mat)
    logger.info("%s raw interaction contact matrix saved to %s_aggLoops.npz" %
                (mat.shape[0], output))
    #1D signal
    if oneD:
        ss = pd.concat([d[1] for d in ds if d[1] is not None], axis=0)
        total = meta["Unique PETs"] * 2
        #normalized to RPM
        ss = ss / total * 10**6
        logger.info(
            "%s RPM normalized 1D signal around anchors saved to %s_aggLoops_1D.txt"
            % (ss.shape[0], output))
        ss.to_csv(output + "_aggLoops_1D.txt", sep="\t")
    else:
        ss = None
    nloops = []
    for key in loops.keys():
        nloops.extend(loops[key])
    bwSigs = {}
    if len(bws) > 0:
        logger.info("Getting the 1D bigWig signals.")
        for bw in bws:
            bwName, bwSig = getBwSig(bw, nloops, ext=ext, skipZeros=skipZeros)
            bwSigs[bwName] = bwSig
    else:
        bwSigs = None
    plotALoops(mat, output, ss, bwSigs, norm=norm,vmin=vmin,vmax=vmax)


### aggregated viewpoints releated functions
def getAViewPoints(
        key,
        fixy,
        viewPoints,
        upExt=100000,
        downExt=100000,
        bs=1000,
        cut=0,
        mcut=-1,
        skipZeros=False,
        oneD=False,
        oneBins=200,
):
    """
    Get the interaction scores from the view points. 
    """
    key2, mat = parseIxy(fixy, cut=cut, mcut=mcut)
    if mat.shape[0] == 0:
        print(
            "No PETs found in %s maybe due to distance cutoff for PET > %s <%s."
            % (fixy, cut, mcut))
        return None
    #XY object
    xy = XY(mat[:, 0], mat[:, 1])
    print(
        "Getting the contact matrix for target %s view points and nearby regions from %s."
        % (len(viewPoints), key))
    mats = []
    vps = []
    for vp in tqdm(viewPoints):
        center = (vp[1] + vp[2]) / 2
        start = max(0, center - upExt)
        if start == 0:
            continue
        end = center + downExt
        rs = xy.queryPeakBoth(start, end)
        nmat = mat[list(rs), ]
        nmat = getObsMat(nmat, start, end, bs)
        nmat = nmat[:-1, :-1]
        if skipZeros and np.sum(nmat) == 0:
            continue
        mats.append(nmat)
        #prepare for 1d profiles
        vps.append([vp[0], start, end])
    mats = np.array(mats)
    if oneD:
        ss = get1DSigMat(xy, vps, bins=oneBins)
    else:
        ss = None
    return mats, ss


def getBwSigV(bw,
              viewPoints,
              upExt=100000,
              downExt=100000,
              skipZeros=False,
              bins=200):
    """
    Get the 1D sig for bigwig around view points.
    """
    print("Getting signal from %s around view points." % bw)
    name = bw.split("/")[-1].split(".bw")[0]
    bw = pyBigWig.open(bw)
    ys = np.zeros(bins)
    c = 0
    for vp in tqdm(viewPoints):
        center = int((vp[1] + vp[2]) / 2)
        start = max(0, center - upExt)
        if start == 0:
            continue
        end = center + downExt
        try:
            s = np.array(bw.values(vp[0], start, end))
        except:
            continue
        s = np.nan_to_num(s)
        s = getBinMean(s, bins)
        if skipZeros:
            if np.sum(s) > 0:
                ys = ys + s
                c += 1
        else:
            ys = ys + s
            c += 1
    ys = ys / c
    return name, ys


def plotAViewPoints(
        mat,
        fout,
        upExt,
        downExt,
        bs,
        ss=None,
        bwSigs=None,
        bins=200,
        norm=False,
        vmin=None, 
        vmax=None,
):
    """
    bwSigs: {name:signal array}
    """
    n = mat.shape[0]
    #take care of colormap
    center = None

    #get the enrichment score, defined as viewpoint all interactions / mean( upstream&downstream)
    p = int(upExt / bs)
    es = []
    for i in range(mat.shape[0]):
        nmat = mat[i]
        ts = []
        for j in range(nmat.shape[0]):
            ts.append(nmat[:j, j].sum() + nmat[j, j:].sum())
        if np.sum(ts) == 0:
            continue
        es.append(ts[p] / np.mean(ts))

    mat = np.log2(mat + 1)
    mat = np.mean(mat, axis=0)
    if norm:
        cmap = cmaps["cool"]
        mat = mat / np.mean(mat)
        mat = np.log2(mat)
        center = 0
        label = "normlized by mean"
    else:
        mat = np.log2(mat + 1)
        cmap = cmaps["red"]
        if vmin is None:
            vmin = 0
        center=None
        label = "log2(PETs+1)"

    #1D signal
    if ss is not None:
        ss = np.mean(ss, axis=0)

    title = "%s view points; ES: %.3f" % (n, np.mean(es))

    #prepare fig
    if bwSigs is not None and ss is not None:
        fig = pylab.figure(figsize=(2, 3 + len(bwSigs) * 0.5))
        hr = [0.5] * len(bwSigs)
        hr.extend([0.5, 2.4, 0.1])
        gs = mpl.gridspec.GridSpec(len(bwSigs) + 3, 1, height_ratios=hr)
    elif bwSigs is None and ss is not None:
        fig = pylab.figure(figsize=(2, 3))  #the extra 1, 0.5 for 1D signal,
        hr = [0.5, 2.4, 0.1]
        gs = mpl.gridspec.GridSpec(3, 1, height_ratios=hr)
    elif bwSigs is not None and ss is None:
        fig = pylab.figure(figsize=(2, 2.5 + len(bwSigs) * 0.5))
        hr = [0.5] * len(bwSigs)
        hr.extend([2.4, 0.1])
        gs = mpl.gridspec.GridSpec(len(bwSigs) + 2, 1, height_ratios=hr)
    else:
        fig = pylab.figure(figsize=(2, 2.5))  #the extra 1, 0.5 for 1D signal,
        hr = [2.4, 0.1]
        gs = mpl.gridspec.GridSpec(2, 1, height_ratios=hr)
    fig.suptitle(title, fontsize=8)

    xs = list(range(bins))

    #plot bigWig signals
    if bwSigs is not None:
        for i, name in enumerate(list(bwSigs.keys())):
            ax = fig.add_subplot(gs[i])
            ax.plot(xs, bwSigs[name], color=colors[i], label=name)
            ax.tick_params(axis='both', which='major', labelsize=4)
            ax.set_xticklabels([])
            ax.set_xlim([np.min(xs), np.max(xs)])
            ax.set_ylim([0, np.max(bwSigs[name]) * 1.1])
            ax.legend(fontsize=5, fancybox=False, frameon=False)

    if ss is not None:
        #plot 1D signal
        if bwSigs is not None:
            ax = fig.add_subplot(gs[len(bwSigs)])
            ax.plot(xs, ss, color=colors[len(bwSigs)], label="mean 1D profile")
        else:
            ax = fig.add_subplot(gs[0])
            ax.plot(xs, ss, color=colors[0], label="mean 1D profile")
        ax.legend(fontsize=5, fancybox=False, frameon=False)
        ax.tick_params(axis='both', which='major', labelsize=4)
        ax.set_xlim([np.min(xs), np.max(xs)])
        ax.set_ylim([0, np.max(ss) * 1.1])
        ax.set_xticklabels([])
        ax.set_ylabel("RPM", fontsize=6)

    #plot mean heatmap
    ax = fig.add_subplot(gs[-2])
    cax = fig.add_subplot(gs[-1])
    sns.heatmap(mat,
                xticklabels=False,
                yticklabels=False,
                square=True,
                ax=ax,
                cmap=cmap,
                cbar_ax=cax,
                linewidths=0,
                vmin=vmin,
                vmax=vmax,
                center=center,
                cbar_kws={
                    'orientation': 'horizontal',
                    "shrink": 0.2,
                })
    ax.text(ax.get_xlim()[0],
            ax.get_ylim()[1] * 1.2,
            "-%.2f kb" % (upExt / 1000),
            fontsize=5)
    ax.text(ax.get_xlim()[1] * 0.8,
            ax.get_ylim()[1] * 1.2,
            "%.2f kb" % (downExt / 1000),
            fontsize=5)
    ax.text(p * 0.8, ax.get_ylim()[1] * 1.2, "view points", fontsize=5)
    cax.tick_params(labelsize=4)
    ax.axvline(x=ax.get_xlim()[0], color="k", linewidth=2)
    ax.axvline(x=ax.get_xlim()[1], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[0], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[1], color="k", linewidth=2)
    pylab.savefig("%s_aggViewPoints.pdf" % fout)


def aggViewPoints(
        predir,
        viewPointF,
        output,
        logger,
        bws="",
        upExt=100000,
        downExt=100000,
        bs=5000,
        cut=0,
        mcut=-1,
        cpu=1,
        skipZeros=False,
        oneD=False,
        oneBins=500,
        norm=False,
        vmin=None,
        vmax=None,
):
    """
    Aggregated view points analysis. 
    """
    if bws != "":
        bws = bws.split(",")
        bws = [bw for bw in bws if os.path.isfile(bw)]
    else:
        bws = []
    #read in view points orgainzed by chromosomes
    viewPoints = {}
    for line in open(viewPointF):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        key = line[0] + "-" + line[0]
        if key not in viewPoints:
            viewPoints[key] = []
        line[1] = int(line[1])
        line[2] = int(line[2])
        viewPoints[key].append(line)
    #meta data
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(viewPoints.keys())))
    #get the data
    ds = Parallel(n_jobs=cpu,
                  backend="multiprocessing")(delayed(getAViewPoints)(
                      key,
                      meta["data"]["cis"][key]["ixy"],
                      viewPoints[key],
                      upExt=upExt,
                      downExt=downExt,
                      bs=bs,
                      cut=cut,
                      mcut=mcut,
                      skipZeros=skipZeros,
                      oneD=oneD,
                      oneBins=oneBins,
                  ) for key in keys)
    #re-complie results
    #matrix
    mat = np.concatenate([d[0] for d in ds if d[0] is not None], axis=0)
    np.savez(output + "_aggViewPoints.npz", mat)
    logger.info(
        "%s raw interaction contact matrix saved to %s_aggViewPoints.npz" %
        (mat.shape[0], output))

    #1D signal
    if oneD:
        ss = pd.concat([d[1] for d in ds if d[1] is not None], axis=0)
        total = meta["Unique PETs"] * 2
        #normalized to RPM
        ss = ss / total * 10**6
        logger.info(
            "%s RPM normalized 1D signal view points saved to %s_aggViewPoints_1D.txt"
            % (ss.shape[0], output))
        ss.to_csv(output + "_aggViewPoints_1D.txt", sep="\t")
    else:
        ss = None

    if len(bws) > 0:
        vps = []
        for k, v in viewPoints.items():
            vps.extend(v)
        logger.info("Getting the 1D bigWig signals.")
        bwSigs = {}
        for bw in bws:
            bwName, bwSig = getBwSigV(
                bw,
                vps,
                upExt=upExt,
                downExt=downExt,
                skipZeros=skipZeros,
                bins=oneBins,
            )
            bwSigs[bwName] = bwSig
    else:
        bwSigs = None

    #plot the data
    plotAViewPoints(mat,
                    output,
                    upExt,
                    downExt,
                    bs,
                    ss,
                    bwSigs,
                    bins=oneBins,
                    norm=norm,
                    vmin=vmin,
                    vmax=vmax,
                    )


### aggregated two anchors releated functions
def getATwoAnchors(key,
                   fixy,
                   loops,
                   ext=0.25,
                   bs=100,
                   cut=0,
                   mcut=-1,
                   skipZeros=False,
                   oneBins=200,
                   oneD=False):
    """
    Get the contact matrix for target regions.
    """
    key2, mat = parseIxy(fixy, cut=cut, mcut=mcut)
    if mat.shape[0] == 0:
        print(
            "No PETs found in %s maybe due to distance cutoff for PET > %s <%s."
            % (fixy, cut, mcut))
        return None, None
    xy = XY(mat[:, 0], mat[:, 1])
    print("Get the contact matrix nearby target %s loops from %s." %
          (len(loops), key))
    mats = []
    nrs = []
    nloops = []
    cabs = []
    for loop in tqdm(loops):
        start = min(loop.x_start, loop.y_start)
        end = max(loop.x_end, loop.y_end)
        size = end - start
        start = start - ext * size
        end = end + ext * size
        if start < 0:
            continue
        bbs = int((end - start) / bs)
        mat = np.zeros([bs, bs])
        rs = []
        for i in range(bs):
            rs.append([start + bbs * i, start + bbs * (i + 1)])
        ns = []
        for iv in rs:
            ns.append(xy.queryPeak(iv[0], iv[1]))
        for i, na in enumerate(ns):
            for j, nb in enumerate(ns):
                nrab = float(len(na.intersection(nb)))
                mat[i][j] = nrab
        if skipZeros and np.sum(mat) == 0:
            continue
        ca, cb, cab = xy.queryLoop(loop.x_start, loop.x_end, loop.y_start,
                                   loop.y_end)
        cab = float(
            len(cab)) / (loop.x_end - loop.x_start + loop.y_end - loop.y_start)
        cabs.append(cab)
        nloops.append(loop)
        mats.append(mat)
        nrs.append([loop.chromX, int(start), int(end)])
    mats = np.array(mats)
    if oneD:
        ss = get1DSigMat(xy, nrs, bins=oneBins)
        #reverse some of them
        for i in range(ss.shape[0]):
            if nloops[i].x_center > nloops[i].y_center:
                s = ss[i]
                s = s[::-1]
                ss[i] = s
    else:
        ss = None
    return mats, ss, cabs


def getBwSigT(
        bw,
        loops,
        ext=0.25,
        skipZeros=False,
        bins=200,
):
    """
    Get the 1D sig for bigwig around view points.
    """
    print("Getting signal from %s around view points." % bw)
    name = bw.split("/")[-1].split(".bw")[0]
    bw = pyBigWig.open(bw)
    ys = np.zeros(bins)
    c = 0
    for loop in tqdm(loops):
        start = min(loop.x_start, loop.y_start)
        end = max(loop.x_end, loop.y_end)
        size = end - start
        center = (start + end) / 2
        d = (end - start) * (1 + ext * 2) / 2
        start = int(center - d)
        end = int(center + d)
        #start = int(start - ext * size)
        #end = int(end + ext * size)
        try:
            s = np.array(bw.values(loop.chromX, start, end))
        except:
            continue
        s = np.nan_to_num(s)
        s = getBinMean(s, bins)
        if loop.x_center > loop.y_center:
            s = s[::-1]
        if skipZeros:
            if np.sum(s) != 0:
                ys = ys + s
                c += 1
        else:
            ys = ys + s
            c += 1
    ys = ys / c
    return name, ys


def plotATwoAnchors(mat,
                    cabs,
                    fout,
                    ss=None,
                    bwSigs=None,
                    vmin=None,
                    vmax=None,
                    bins=200):
    """
    Plot the mean border heatmap.
    """
    #1D signal
    if ss is not None:
        ss = np.mean(ss, axis=0)

    mat = np.mean(mat, axis=0)
    mat = np.log2(mat + 1)
    cmap = sns.light_palette("red", n_colors=9).as_hex()
    cmap[0] = "#FFFFFF"
    cmap = ListedColormap(cmap)
    #cmap = cmaps["red"]
    if vmin is None:
        vmin = 0
    center = None
    label = "log2(PETs+1)"

    title = "%s loops; %.2f RPKM" % (len(cabs), np.mean(cabs))

    #prepare fig
    if bwSigs is not None and ss is not None:
        fig = pylab.figure(figsize=(2, 2.8 + len(bwSigs) * 0.3))
        hr = [0.3] * len(bwSigs)
        hr.extend([0.3, 2.4, 0.1])
        gs = mpl.gridspec.GridSpec(len(bwSigs) + 3, 1, height_ratios=hr)
    elif bwSigs is None and ss is not None:
        fig = pylab.figure(figsize=(2, 2.8))  #the extra 1, 0.5 for 1D signal,
        hr = [0.3, 2.4, 0.1]
        gs = mpl.gridspec.GridSpec(3, 1, height_ratios=hr)
    elif bwSigs is not None and ss is None:
        fig = pylab.figure(figsize=(2, 2.5 + len(bwSigs) * 0.3))
        hr = [0.3] * len(bwSigs)
        hr.extend([2.4, 0.1])
        gs = mpl.gridspec.GridSpec(len(bwSigs) + 2, 1, height_ratios=hr)
    else:
        fig = pylab.figure(figsize=(2, 2.5))  #the extra 1, 0.5 for 1D signal,
        hr = [2.4, 0.1]
        gs = mpl.gridspec.GridSpec(2, 1, height_ratios=hr)
    #gs.update(hspace=0.1)
    fig.suptitle(title, fontsize=8)

    xs = list(range(bins))

    #plot bigWig signals
    if bwSigs is not None:
        for i, name in enumerate(list(bwSigs.keys())):
            ax = fig.add_subplot(gs[i])
            ax.plot(xs, bwSigs[name], color=colors[i], label=name)
            ax.tick_params(axis='both', which='major', labelsize=4)
            ax.set_xticklabels([])
            ax.set_xlim([np.min(xs), np.max(xs)])
            ax.set_ylim([0, np.max(bwSigs[name]) * 1.1])
            ax.legend(fontsize=5, fancybox=False, frameon=False)

    if ss is not None:
        #plot 1D signal
        if bwSigs is not None:
            ax = fig.add_subplot(gs[len(bwSigs)])
            ax.plot(xs, ss, color=colors[len(bwSigs)], label="mean 1D profile")
        else:
            ax = fig.add_subplot(gs[0])
            ax.plot(xs, ss, color=colors[0], label="mean 1D profile")
        ax.legend(fontsize=5, fancybox=False, frameon=False)
        ax.tick_params(axis='both', which='major', labelsize=4)
        ax.set_xlim([np.min(xs), np.max(xs)])
        ax.set_ylim([0, np.max(ss) * 1.1])
        ax.set_xticklabels([])
        ax.set_ylabel("RPM", fontsize=6)

    #plot mean heatmap
    ax = fig.add_subplot(gs[-2])
    cax = fig.add_subplot(gs[-1])
    sns.heatmap(mat,
                xticklabels=False,
                yticklabels=False,
                square=True,
                ax=ax,
                cmap=cmap,
                cbar_ax=cax,
                linewidths=0,
                vmin=vmin,
                center=center,
                vmax=vmax,
                cbar_kws={
                    'orientation': 'horizontal',
                    "shrink": 0.2,
                })
    cax.tick_params(labelsize=4)
    ax.axvline(x=ax.get_xlim()[0], color="k", linewidth=2)
    ax.axvline(x=ax.get_xlim()[1], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[0], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[1], color="k", linewidth=2)
    pylab.savefig("%s_aggTwoAnchors.pdf" % fout)


def aggTwoAnchors(
        predir,
        loopf,
        output,
        logger,
        bws="",
        ext=0.25,
        cut=0,
        mcut=-1,
        lcut=0,
        cpu=1,
        skipZeros=False,
        norm=False,
        oneD=False,
        vmin=None,
        vmax=None,
):
    """
    Aggregated two anchors analysis.
    """
    if bws != "":
        bws = bws.split(",")
        bws = [bw for bw in bws if os.path.isfile(bw)]
    else:
        bws = []

    loops = parseTxt2Loops(loopf, cut=lcut)
    if len(loops) == 0:
        logger.info("No loops found in %s with loop ditance > %s, return." %
                    (loopf, lcut))
        return

    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(loops.keys())))
    ds = Parallel(n_jobs=cpu,
                  backend="multiprocessing")(delayed(getATwoAnchors)(
                      key,
                      meta["data"]["cis"][key]["ixy"],
                      loops[key],
                      ext=ext,
                      cut=cut,
                      mcut=mcut,
                      skipZeros=skipZeros,
                      oneD=oneD,
                  ) for key in keys)

    #matrix
    mat = np.concatenate([d[0] for d in ds if d[0] is not None], axis=0)
    np.savez(output + "_aggTwoAnchors.npz", mat)
    logger.info(
        "%s raw interaction contact matrix saved to %s_aggTwoAnchors.npz" %
        (mat.shape[0], output))

    #1D signal
    if oneD:
        ss = pd.concat([d[1] for d in ds if d[1] is not None], axis=0)
        total = meta["Unique PETs"] * 2
        #normalized to RPM
        ss = ss / total * 10**6
        logger.info(
            "%s RPM normalized 1D signal around anchors saved to %s_aggTwoAnchors_1D.txt"
            % (ss.shape[0], output))
        ss.to_csv(output + "_aggTwoAnchors_1D.txt", sep="\t")
    else:
        ss = None
    #loop densities
    cabs = np.concatenate([d[2] for d in ds if d[2] is not None], axis=0)
    total = meta["Unique PETs"]
    cabs = cabs / total * 10**9

    nloops = []
    for key in loops.keys():
        nloops.extend(loops[key])
    bwSigs = {}
    if len(bws) > 0:
        logger.info("Getting the 1D bigWig signals.")
        for bw in bws:
            bwName, bwSig = getBwSigT(bw, nloops, ext=ext, skipZeros=skipZeros)
            bwSigs[bwName] = bwSig
    else:
        bwSigs = None
    plotATwoAnchors(mat,
                    cabs,
                    output,
                    ss=ss,
                    bwSigs=bwSigs,
                    vmin=vmin,
                    vmax=vmax)


### aggregated domains releated funcitons
def getDomainBwSig(bw, domains, bins=200, ext=0.5, skipZeros=False):
    """
    Get the 1D sig for bigwig.
    """
    print("Getting signal from %s around domains." % bw)
    name = bw.split("/")[-1].split(".bw")[0]
    bw = pyBigWig.open(bw)
    ys = np.zeros(bins)
    c = 0
    for domain in tqdm(domains):
        size = domain[2] - domain[1]
        start = domain[1] - ext * size
        end = domain[2] + ext * size
        start = int(start)
        end = int(end)
        try:
            sx = bw.values(domain[0], start, end)
        except:
            continue
        sx = np.nan_to_num(sx)
        sx = getBinMean(sx, bins)
        if skipZeros and np.sum(sx) == 0.0:
            continue
        ys = ys + sx
        c += 1
    ys = ys / c
    return name, ys


def getADomains(key,
                fixy,
                domains,
                ext=0.5,
                bs=100,
                cut=0,
                mcut=-1,
                skipZeros=False,
                oneBins=200,
                oneD=False):
    """
    Get the contact matrix for target regions.
    """
    key2, mat = parseIxy(fixy, cut=cut, mcut=mcut)
    if mat.shape[0] == 0:
        print(
            "No PETs found in %s maybe due to distance cutoff for PET > %s <%s."
            % (fixy, cut, mcut))
        return None, None
    xy = XY(mat[:, 0], mat[:, 1])
    print("Get the contact matrix nearby target %s domains from %s." %
          (len(domains), key))
    mats = []
    nrs = []
    es = []
    for domain in tqdm(domains):
        t = xy.queryPeak(domain[1], domain[2])
        b = xy.queryPeakBoth(domain[1], domain[2])
        n = t.difference(b)
        if len(n) > 0:
            e = len(b) / float(len(n))
            es.append(e)
        size = domain[2] - domain[1]
        start = domain[1] - ext * size
        end = domain[2] + ext * size
        if start < 0:
            continue
        bbs = int((end - start) / bs)
        mat = np.zeros([bs, bs])
        rs = []
        for i in range(bs):
            rs.append([start + bbs * i, start + bbs * (i + 1)])
        ns = []
        for iv in rs:
            ns.append(xy.queryPeak(iv[0], iv[1]))
        for i, na in enumerate(ns):
            for j, nb in enumerate(ns):
                nrab = float(len(na.intersection(nb)))
                mat[i][j] = nrab
        if skipZeros and np.sum(mat) == 0:
            continue
        mats.append(mat)
        nrs.append([domain[0], int(start), int(end)])
    mats = np.array(mats)
    if oneD:
        ss = get1DSigMat(xy, nrs, bins=oneBins)
    else:
        ss = None
    return mats, ss, es


def plotADomains(mat,
                 fout,
                 es,
                 size,
                 ss=None,
                 bwSigs=None,
                 vmin=None,
                 vmax=None):
    """
    Plot the mean border heatmap.
    """
    es = np.mean(es)

    #1D signal
    if ss is not None:
        ss = np.mean(ss, axis=0)

    #simple normalization
    n = mat.shape[0]
    label = "log2(PETs+1)"
    mat = np.mean(mat, axis=0)
    mat = np.log2(mat + 1)

    #take the uppper matrix
    mat = rotate(mat, angle=45)
    to = int(mat.shape[0] / 2) - 1
    mat = mat[:to, ]

    #take care of colormap
    #cmap = sns.cubehelix_palette(light=1, as_cmap=True)
    #following ligh color much better
    cmap = sns.light_palette("red", n_colors=9).as_hex()
    cmap[0] = "#FFFFFF"
    cmap = ListedColormap(cmap)
    if vmin == None:
        vmin = 0

    title = "%s domains;ES:%.2f;median size:%.2fkb" % (n, es,
                                                       np.median(size) / 1000)
    #only domain heatmap
    if bwSigs is None and ss is None:
        fig, ax = pylab.subplots(figsize=(3, 1.5))
        sns.heatmap(mat,
                    xticklabels=False,
                    yticklabels=False,
                    square=False,
                    linewidths=0.0,
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={
                        'label': label,
                        "shrink": 0.5
                    })

        ax.set_title(title)
    else:
        #prepare fig
        if bwSigs is not None and ss is not None:
            figsize = (3, 1.5 + 0.5 + len(bwSigs) * 0.5)
            hr = [1.5, 0.5]
            hr.extend([0.5] * len(bwSigs))
        elif bwSigs is None and ss is not None:
            figsize = (3, 2)
            hr = [1.5, 0.5]
        elif bwSigs is not None and ss is None:
            figsize = (3, 1.5 + len(bwSigs) * 0.5)
            hr = [1.5]
            hr.extend([0.5] * len(bwSigs))

        fig, axs = pylab.subplots(len(hr),
                                  1,
                                  gridspec_kw={'height_ratios': hr},
                                  figsize=figsize)

        fig.suptitle(title)
        #heatmap
        ax = axs[0]
        cax = inset_axes(ax, width="3%", height="40%", loc="right")
        sns.heatmap(mat,
                    xticklabels=False,
                    yticklabels=False,
                    square=False,
                    linewidths=0.0,
                    cmap=cmap,
                    ax=ax,
                    cbar_ax=cax,
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={
                        'label': label,
                    })
        cax.tick_params(labelsize=4)
        #plot bigWig signals
        if bwSigs is not None:
            for i, name in enumerate(list(bwSigs.keys())):
                ax = axs[i + 1]
                xs = list(range(len(bwSigs[name])))
                ax.plot(xs, bwSigs[name], color=colors[i], label=name)
                ax.tick_params(axis='both', which='major', labelsize=4)
                ax.set_xticklabels([])
                ax.set_xlim([np.min(xs), np.max(xs)])
                #ax.set_ylim([0, np.max(bwSigs[name]) * 1.1])
                ax.set_ylim(
                    [np.min(bwSigs[name]) * 0.9,
                     np.max(bwSigs[name]) * 1.1])
                ax.legend(fontsize=5, fancybox=False, frameon=False)
        #plot oneD data
        if ss is not None:
            if bwSigs is not None:
                ind = 1 + len(bwSigs)
            else:
                ind = 1
            ax = axs[ind]
            xs = np.arange(len(ss))
            ax.plot(xs, ss, color=colors[ind], label="mean 1D profile")
            ax.legend(fontsize=5, fancybox=False, frameon=False)
            ax.tick_params(axis='both', which='major', labelsize=4)
            ax.set_xlim([np.min(xs), np.max(xs)])
            ax.set_ylim([0, np.max(ss) * 1.1])
            ax.set_xticklabels([])
            ax.set_ylabel("RPM", fontsize=6)

    pylab.savefig("%s_aggDomains.pdf" % fout)


def aggDomains(predir,
               domainf,
               output,
               logger,
               bws="",
               ext=0.5,
               cut=0,
               mcut=-1,
               cpu=1,
               vmin=None,
               vmax=None,
               skipZeros=False,
               oneD=False):
    """
    Aggregated domains analysis. 
    """
    if bws != "":
        bws = bws.split(",")
        bws = [bw for bw in bws if os.path.isfile(bw)]
    else:
        bws = []

    #domains
    size = []
    domains = {}
    for line in open(domainf):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        key = line[0] + "-" + line[0]
        if key not in domains:
            domains[key] = []
        line[1] = int(line[1])
        line[2] = int(line[2])
        domains[key].append(line)
        size.append(line[2] - line[1])

    #pre data
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(domains.keys())))

    #get
    ds = Parallel(n_jobs=cpu, backend="multiprocessing")(delayed(getADomains)(
        key,
        meta["data"]["cis"][key]["ixy"],
        domains[key],
        ext=ext,
        cut=cut,
        mcut=mcut,
        oneD=oneD,
    ) for key in keys)

    mat = np.concatenate([d[0] for d in ds if d[0] is not None], axis=0)
    np.savez(output + "_aggDomains.npz", mat)
    logger.info(
        "%s raw interaction contact matrix saved to %s_aggDomains.npz" %
        (mat.shape[0], output))

    #1D signal
    if oneD:
        ss = pd.concat([d[1] for d in ds if d[1] is not None], axis=0)
        total = meta["Unique PETs"] * 2
        #normalized to RPM
        ss = ss / total * 10**6
        logger.info(
            "%s RPM normalized 1D signal around domains saved to %s_aggDomains_1D.txt"
            % (ss.shape[0], output))
        ss.to_csv(output + "_aggDomains_1D.txt", sep="\t")
    else:
        ss = None

    #enrichment score
    es = np.concatenate([d[2] for d in ds])

    #get bigWig signals
    ndomains = []
    for key in domains.keys():
        ndomains.extend(domains[key])
    bwSigs = {}
    if len(bws) > 0:
        print("Getting the 1D bigWig signals.")
        for bw in bws:
            bwName, bwSig = getDomainBwSig(bw, ndomains, ext=ext)
            bwSigs[bwName] = bwSig
    else:
        bwSigs = None

    #plot the aggregated domains
    plotADomains(mat, output, es, size, ss, bwSigs, vmin, vmax)
