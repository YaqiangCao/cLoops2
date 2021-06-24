#!/usr/bin/env python
#--coding:utf-8 --
"""
plot.py
Plotting related functions for cLoops2.
2020-02-25: refine plot bigWig tracks, reduce the figure size by bin the data
2020-06-24: going to add virtual 4C plot such as https://www.researchgate.net/publication/317638806_Evolutionary_Analysis_of_Candidate_Non-Coding_Elements_Regulating_Neurodevelopmental_Genes_in_Vertebrates/figures?lo=1 
2020-08-25: adding eigenvector for the plot of heatmap, not available in command line
2020-08-27: adding genes for the plot of heatmap
2020-09-03: change line style of bed annotations to recetangle style.
2020-09-13: refine virtual4C plot
2020-10-30: adding plotting PETs as arches
2021-02-04: adding option for not requring heatmps/arches
2021-02-09: refine some details of code
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os
import json
import argparse
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import sparse
import matplotlib as plt
from matplotlib.patches import Arc
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from sklearn.decomposition import PCA
from matplotlib.ticker import AutoLocator
from matplotlib.colors import ListedColormap
from scipy.ndimage import rotate  #rotate the domain

#cLoops2
from cLoops2.ds import XY, Exon, Gene
from cLoops2.settings import *
from cLoops2.utils import getLogger
from cLoops2.io import parseIxy, parseTxt2Loops
from cLoops2.cmat import getObsMat, getExpMat, get1DSig, getBinMean, getVirtual4CSig


def plotGmmEst(dis, ps, eps, fout):
    """
    Plot the distribution of distances estimated from the gaussian mixture models. 
    @param dis: numpy.array, distance vector
    @param ps: numpy.array, classification of points type vector
    @param eps: int, eps
    @param fout: str, outpuf pdf file name
    """
    fig, axs = pylab.subplots(1, 2, figsize=(4, 2.75), sharey=True)
    #raw PETs
    sns.kdeplot(dis, ax=axs[0], shade=True, color=colors[2])
    axs[0].set_ylabel("Density")
    axs[0].set_title("Raw")
    #gmm infered PETs
    nsa = np.where(ps == 0)[0]
    nsb = np.where(ps == 1)[0]
    if np.mean(dis[nsa]) > np.mean(
            dis[nsb]):  #change the order if first classes mean larger
        nsa, nsb = nsb, nsa
    sns.kdeplot(dis[nsa],
                ax=axs[1],
                shade=True,
                color=colors[0],
                label="potential peak PETs")  #colors from the cLoops2.settings
    sns.kdeplot(dis[nsb],
                ax=axs[1],
                shade=True,
                color=colors[1],
                label="potential loop PETs")  #colors from the cLoops2.settings
    axs[1].legend()
    axs[1].set_title("GMM inferred PETs\nEstimated eps=%s" % eps)
    #set common x-label
    fig.text(0.5, 0.0, "Distance between two ends (log2,bp)", ha='center')
    pylab.tight_layout()
    pylab.savefig(fout)


def plotKDis(dis, k, fout):
    """
    Plot the k-distance distribution. 
    @param dis: numpy.array, distance vector
    @param k: int, k-neighbor
    @param fout: str, outpuf pdf file name
    """
    fig, ax = pylab.subplots()
    x = np.arange(len(dis))
    #ax.scatter( x,dis,color=colors[0],s=1 )
    ax.plot(x, dis, color=colors[0])
    ax.set_xlabel("Points sorted by distance")
    #ax.set_ylabel("%s-NN distance (log2,bp)" % k)
    ax.set_ylabel("%s-NN distance (log10,bp)" % k)
    #ax.set_xscale("log")
    pylab.tight_layout()
    pylab.savefig(fout)


def plotKDisE(dis, k, knee, eps, fout):
    """
    Plot the k-distance distribution, enhanced, if can auto detect the knee.
    @param dis: numpy.array, distance vector
    @param k: int, k-neighbor
    @param knee:int, the knee for the plot, show the vertical line 
    @param eps: float, the log2 transform eps, show the horizontal 
    @param fout: str, outpuf pdf file name
    """
    fig, ax = pylab.subplots()
    x = np.arange(len(dis))
    #ax.scatter( x,dis,color=colors[0],s=1 )
    ax.plot(x, dis, color=colors[0])
    ax.set_xlabel("Points sorted by distance")
    #ax.set_ylabel("%s-NN distance (log2,bp)" % k)
    ax.set_ylabel("%s-NN distance" % k)
    ax.axvline(knee, label="knee", linestyle="--", color=colors[1])
    ax.axhline(eps,
               label="estimated eps:%s bp" % (int(2**eps)),
               linestyle="--",
               color=colors[2])
    ax.legend()
    #ax.set_xscale("log")
    ax.set_xticklabels([])
    pylab.tight_layout()
    pylab.savefig(fout)


def plotIntraCut(di, ds, cut, log=True, prefix="test"):
    """
    Plot the distance cutoff of self-ligation and inter-ligation reads.
    """
    di = np.abs(np.array(di))
    ds = np.abs(np.array(ds))
    di = di[~np.isnan(di)]
    ds = ds[~np.isnan(ds)]
    di = di[di > 0]
    ds = ds[ds > 0]
    if log:
        di = np.log2(di)
        ds = np.log2(ds)
    fig, ax = pylab.subplots()
    sns.kdeplot(di,
                ax=ax,
                shade=True,
                label="inter-ligation PETs:%s" % len(di),
                color=colors[0])
    sns.kdeplot(ds,
                ax=ax,
                shade=True,
                label="self-ligation PETs:%s" % len(ds),
                color=colors[1])
    ax.axvline(np.log2(cut),
               label="distance cutoff:%.2f kb" % (cut / 1000.0),
               color=colors[2])
    leg = ax.legend(loc="best", shadow=True, fancybox=True)
    ax.set_xlabel("Distance between PET ends (log2(bp))")
    ax.set_ylabel("Density")
    pylab.savefig("%s.pdf" % prefix)


def plotEstRes(binSizes, cumBins, singletonRatios, PETsRatio, prefix):
    """
    Plot the estimation of resolution.
    """
    fig, ax = pylab.subplots()
    for i in range(len(binSizes)):
        ss = cumBins[i]
        r = singletonRatios[i]
        tp = int(round(r))
        pr = PETsRatio[i]
        ax.plot(ss.index[:tp],
                ss.values[:tp],
                color=colors[i],
                linestyle="--",
                linewidth=1)
        ax.plot(ss.index[tp:],
                ss.values[tp:],
                label="%s: %.2f %% bins %.2f %% PETs" %
                (binSizes[i], 100 - r, 100 - pr),
                color=colors[i],
                linewidth=2)
    ax.plot([0, 100], [0, 100], color="k", label="random")
    ax.set_xlim([0, 100])
    ax.set_ylim([0, 100])
    ax.legend(fontsize=6, fancybox=False, frameon=False)
    ax.set_xlabel("Percentage of contact matrix bins")
    ax.set_ylabel("Percentage of PETs")
    pylab.savefig("%s_estRes.pdf" % prefix)


def plotEstSat(binSizes, totPETs, data, tol, prefix):
    """
    Plot the estimation of sequencing signal saturation.
    """
    fig, ax = pylab.subplots()
    for i in range(len(binSizes)):
        d = data[i]
        xs = d.index
        ys = d.mean(axis=1)
        std = d.std(axis=1)
        ax.errorbar(xs,
                    ys,
                    yerr=std,
                    capsize=3,
                    linewidth=1,
                    elinewidth=1,
                    color=colors[i],
                    label="resolution:%s" % binSizes[i])
        #ax.set_xticks(xs)
    ax.legend()
    ax.set_title("total PETs:%s M" % (totPETs / 10**6))
    ax.set_xlabel("Sub-sampling ratio")
    ax.set_ylabel("Detected contact matrix bins ratio (>=%sPETs)" % tol)
    pylab.savefig("%s_estSat.pdf" % prefix)


def plotCorrScatterPCC(mat, fout):
    """
    Density scatter plot for two samples correlation
    """
    fig, ax = pylab.subplots()
    cmap = sns.light_palette("red", n_colors=9).as_hex()
    cmap[0] = "#FFFFFF"
    cmap = ListedColormap(cmap)
    #cmap = sns.cubehelix_palette(light=1, as_cmap=True)
    #print(type(cmap))
    da = mat.columns[0]
    db = mat.columns[1]
    sa = mat[da]
    sb = mat[db]
    corr = sa.corr(sb)
    hb = ax.hexbin(sa + 1,
                   sb + 1,
                   gridsize=100,
                   cmap=cmap,
                   bins="log",
                   xscale="log",
                   yscale="log")
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('log10(N), number of points')
    ax.set_title("vector size:%s\nPCC:%.3f" % (len(sa), corr))
    ax.set_xlabel(da + ",log10(PET+1)")
    ax.set_ylabel(db + ",log10(PET+1)")
    pylab.savefig("%s" % fout)


def plotCorrScatterPCA(mat, fout):
    """
    Density scatter plot for two samples correlation
    """
    fig, ax = pylab.subplots()
    cmap = sns.light_palette("red", n_colors=9).as_hex()
    cmap[0] = "#FFFFFF"
    cmap = ListedColormap(cmap)
    da = mat.columns[0]
    db = mat.columns[1]
    sa = mat[da]
    sb = mat[db]
    corr = sa.corr(sb)
    hb = ax.hexbin(
        sa,
        sb,
        gridsize=100,
        cmap=cmap,
        bins="log",
    )
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('log10(N), number of points')
    ax.set_title("vector size:%s\nPCC:%.3f" % (len(sa), corr))
    ax.set_xlabel(da + ", top PCs")
    ax.set_ylabel(db + ", top PCs")
    pylab.savefig("%s" % fout)


def plotCorrHeatmap(mat, fout):
    """
    Correlation heatmap plot for two samples correlation.
    """
    #fig, ax = pylab.subplots(
    cmap = sns.diverging_palette(250, 15, s=75, l=40, n=11).as_hex()
    cmap[int(len(cmap) / 2)] = "#FFFFFF"
    cmap = ListedColormap(cmap)
    g = sns.clustermap(
        mat,
        xticklabels=False,
        yticklabels=True,
        square=True,
        center=0,
        linewidths=0.0,
        cmap=cmap,
        figsize=(0.5 * mat.shape[1], 0.5 * mat.shape[1]),
        annot=True,
        fmt=".3f",
        annot_kws={
            "size": "3",
            'label': "PCC",
        },
    )
    pylab.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    pylab.savefig(fout)


def getBedRegion(f, chrom, start, end):
    """
    Get the target region in bed file.
    """
    rs = []
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        try:
            c = line[0]
            s = int(line[1])
            e = int(line[2])
        except:
            cotninue
        if c != chrom:
            continue
        if s >= start and e <= end:
            rs.append([s, e])
    return rs


def parseGtf(line):
    """
    Parse gene gtf line.
    """
    e = Exon()
    e.chrom = line[0]
    e.start = int(line[3])
    e.end = int(line[4])
    e.length = e.end - e.start
    e.strand = line[6]
    attr = line[8].replace('"', '').split(";")
    ts = {}
    for t in attr:
        t = t.split()
        if len(t) != 2:
            continue
        ts[t[0]] = t[1]
    e.name = ts["gene_name"]
    e.id = ts["gene_id"]
    return e


def stichExons(exons, margin=1):
    """
    Stich close exons based on postion array. 
    """
    cov = set()
    for i, exon in enumerate(exons):
        cov.update(range(exon.start, exon.end + 1))
    cov = list(cov)
    cov.sort()
    nexons = []
    i = 0
    while i < len(cov) - 1:
        for j in range(i + 1, len(cov)):
            if cov[j] - cov[j - 1] > margin:
                break
            else:
                continue
        exon = exons[0].chrom
        exon = Exon()
        exon.chrom = exons[0].chrom
        exon.start = cov[i]
        exon.end = cov[j - 1]
        exon.strand = exons[0].strand
        exon.length = cov[j - 1] - cov[i] + 1
        exon.id = exons[0].id
        exon.name = exons[0].name
        nexons.append(exon)
        i = j  #update search start
    return nexons


def getGenes(f, chrom, start, end):
    """
    Get the target gene in the gtf file.
    """
    gs = {}
    for line in open(f):
        if line.startswith("#"):
            continue
        line = line.split("\n")[0].split("\t")
        if line[0] != chrom:
            continue
        if line[2] != "exon":
            continue
        e = parseGtf(line)
        if e.name not in gs:
            g = Gene()
            g.chrom = e.chrom
            g.start = e.start
            g.end = e.end
            g.strand = e.strand
            g.name = e.name
            g.id = e.id
            g.exons = {(e.start, e.end): e}
            gs[g.name] = g
        else:
            #same position exons
            if (e.start, e.end) in gs[e.name].exons:
                continue
            else:
                g = gs[e.name]
                if e.start < g.start:
                    g.start = e.start
                if e.end > g.end:
                    g.end = e.end
                g.exons[(e.start, e.end)] = e
    #select genes in the target region
    ngs = {}
    for n, g in gs.items():
        if g.start < start or g.end > end:
            continue
        g.exons = stichExons(list(g.exons.values()))
        ngs[n] = g
    return ngs


def parseBwvs(bws, bwvs=""):
    """
    Parse input bigwig values limts.
    """
    if bwvs == "":
        bwvs = []
        for i in range(len(bws)):
            bwvs.append([None, None])
    else:
        bwvs = bwvs.split(";")
        nbwvs = []
        for t in bwvs:
            if t == "":
                nbwvs.append([None, None])
            else:
                t = t.split(",")
                t = list(map(float, t))
                t.sort()
                nbwvs.append(t)
        bwvs = nbwvs
    return bwvs


def plotGene(ax, n, g, start, end, space=0.02):
    """
    Plot one genes.
    @param ax: maplotlib ax
    @param n: str, gene name
    @param g: cLoops2:ds:Gene object
    @param start: int, start region for plotting
    @param end: int, end region for plotting
    @param space: float, name and gene distance releative
    """
    ax.axis("off")
    ax.set_xlim([start, end])
    ax.set_ylim([0, 1])
    #plot intron as line, exon as block
    for i, exon in enumerate(g.exons):
        c = "k"
        if g.strand == "+" and i == 0:
            c = colors[1]
        if g.strand == "-" and i == len(g.exons) - 1:
            c = colors[3]
        p = patches.Rectangle((exon.start, 0.1),
                              exon.end - exon.start,
                              0.8,
                              fill=True,
                              color=c,
                              alpha=1)
        ax.add_patch(p)
        if i > 0:
            ax.plot([g.exons[i - 1].end, exon.start], [0.5, 0.5],
                    color="gray",
                    linewidth=0.5,
                    linestyle="--")
    #plot direction and name
    if len(g.exons) > 1:
        if g.strand == "+":
            #c = "green"
            c = colors[1]
            ax.plot([g.exons[0].end, g.exons[1].start], [0.5, 0.5],
                    color=c,
                    linewidth=1,
                    linestyle="-")
            p = g.exons[0].start - (end - start) * (space * 2)
            ax.text(p, 0.15, n, color=c, fontsize=5)
        else:
            #c = "red"
            c = colors[3]
            ax.plot([g.exons[-2].end, g.exons[-1].start], [0.5, 0.5],
                    color=c,
                    linewidth=1,
                    linestyle="-")
            p = g.exons[-1].end + (end - start) * space
            ax.text(p, 0.15, n, color=c, fontsize=5)
    else:
        if g.strand == "+":
            c = colors[1]
            p = g.exons[0].start - (end - start) * (space * 2)
            ax.text(p, 0.15, n, color=c, fontsize=5,style="italic")
        else:
            c = colors[3]
            p = g.exons[-1].end + (end - start) * space
            ax.text(p, 0.15, n, color=c, fontsize=5,style="italic")
    return ax


def plotCoverage(ax, ys, colori=1, label="", vmin=None, vmax=None,
                 lencut=1000):
    """
    Plot 1D coverage data.
    @param ax: matplotlib ax
    @param ys: numpy.array, y-axis coverages
    @param colori: int, color index
    @param label: str, name/label for the data
    @param vmin: float, y-axis vmin
    @param vmax: float, y-axis vmax
    @param lencut: int, if the vector of xs/ys is too long, short them by bin averages
    @return ax: matplotlib ax
    """
    if len(ys) > lencut:
        ys = getBinMean(ys, lencut)
    xs = np.arange(len(ys))
    ax.plot(xs, ys, color=colors[colori], label=label, linewidth=0)
    ax.fill_between(np.arange(len(ys)), 0, ys, color=colors[colori], alpha=0.8)
    ax.set_xticklabels([])
    ax.set_xlim([np.min(xs), np.max(xs)])
    #set y-axis lim
    if vmin is None:
        vmin = np.min(ys)
    if vmax is None:
        vmax = np.max(ys)
    vmin = float("%.3f"%vmin)
    vmax = float("%.3f"%vmax)
    p = (vmax - vmin) * 0.15
    ax.set_yticks([vmin, vmax - p])
    ax.set_yticklabels([str(vmin), str(vmax)])
    ax.set_ylim([vmin,vmax])
    ax.tick_params(axis='both', which='major', labelsize=4)
    ax.legend(fontsize=6, fancybox=False, frameon=False)
    return ax


def plotRegion(ax, rs, start, end, colori=1, label=""):
    """
    Plot genomic region.
    @param ax: matplotlib ax
    @param rs: [start,end], both start and end are ints
    @param colori: int, color index
    @param label: str, name/label for the data
    """
    for r in rs:
        p = patches.Rectangle((r[0], 0.2),
                              r[1] - r[0],
                              0.6,
                              fill=True,
                              color=colors[colori],
                              alpha=0.8)
        ax.add_patch(p)
    ax.set_ylim([0, 1])
    ax.text((start + end) / 2, 0.2, label, fontsize=6)
    ax.axis("off")
    ax.set_xlim([start, end])
    return ax


def plotLoops(ax, loops, nchrom, start, end,xy2=None):
    """
    Plot loops as arches
    """
    cabs = []
    nloops = []
    for loop in loops[nchrom]:
        s = min(loop.x_start, loop.y_start)
        e = max(loop.x_end, loop.y_end)
        if start < s and e < end:
            nloops.append(loop)
            #query from the data for number of PETs
            if xy2 is not None:
                ca, cb, cab = xy2.queryLoop(loop.x_start, loop.x_end,
                                        loop.y_start, loop.y_end)
            else:
                cab = [1]
            cabs.append(len(cab))
    #start plot
    ncabs = [c for c in cabs if c > 0]
    if len(ncabs) > 0:
        minCab = np.min(ncabs)
        #modify line width for arches , just in case the line too wide
        lws = [c / minCab for c in cabs]
        if max(lws) > 10:
            lws = [1] * len(lws)
        pa = 0
        pb = 1.0
        ymax = 0
        for i, loop in enumerate(nloops):
            ca = (loop.x_start + loop.x_end) / 2
            cb = (loop.y_start + loop.y_end) / 2
            cc = (ca + cb) / 2
            npa = float(ca - start) / (end - start) * (pb - pa)
            npb = float(cb - start) / (end - start) * (pb - pa)
            npc = float(cc - start) / (end - start) * (pb - pa)
            a = npb - npa  #a is x axis size for eclipse
            b = a / 2  #b is y axis size for eclipse
            if b > ymax:
                ymax = b
            if cabs[i] < 1:
                continue
            ax.add_patch(
                Arc(
                    (npc, 0),
                    a,
                    b,
                    theta1=0,
                    theta2=180,
                    edgecolor=colors[1],
                    lw=lws[i],
                ))
            if xy2 is not None:
                ax.text(npc, b / 2, cabs[i], fontsize=5)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, ymax * 0.6])
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    return ax


def plotMatHeatmap(
        f,
        fo,
        start=0,
        end=-1,
        res=5000,
        cut=0,
        mcut=-1,
        log=False,
        method="obs",
        oneD=False,
        oneDv="",
        corr=False,
        triu=False,
        norm=False,
        bws=[],
        bwvs="",
        bwcs="",
        beds=[],
        loops=None,
        domains="",
        eig=False,
        eig_r=False,
        gtf="",
        virtual4C=False,
        viewStart=-1,
        viewEnd=-1,
        vmin=None,
        vmax=None,
        width=4,
):
    """
    Plot the contact matrix heatmap with 1D tracks or 2D annotations
    """

    chrom, xy = parseIxy(f, cut=cut, mcut=mcut)
    xy2 = XY(xy[:, 0], xy[:, 1])  #XY object

    if start == 0:
        start = np.min(xy)
    if end == -1:
        end = np.max(xy)
    mat = getObsMat(xy, start, end, res)
    bgmat = None
    if method == "obs/exp":
        bgmat = getExpMat(xy, mat.shape, start, end, res)
    if log:
        if bgmat is None:
            ano = "log10(Obs)"
            mat = np.log10(mat + 1)
        else:
            ano = "log10(Obs/Exp)"
            mat = np.log10(mat + 1) - np.log10(bgmat + 1)
    else:
        if bgmat is not None:
            mat = (mat + 1) / (bgmat + 1)
            ano = "Obs/Exp"
        else:
            ano = "Obs"
    if corr:
        ano = ano + " correlation"
        mat = np.corrcoef(mat)
        mat = np.nan_to_num(mat)
    if norm:
        #diag = np.diag(np.diagonal(mat))
        #mat = mat - diag
        m = np.mean(mat)
        s = np.std(mat)
        mat = (mat - m) / s
        ano = ano + " z-socore normalized"
    if oneD:
        predir = os.path.dirname(os.path.realpath(f))
        metaf = predir + "/petMeta.json"
        meta = json.loads(open(metaf).read())
        total = meta["Unique PETs"] * 2
        sig = get1DSig(xy2, start, end)
        sig = sig / total * 10**6
    if eig:
        nmat = np.sum(mat, axis=0).astype("int")
        ps = np.where(nmat == 0)[0]
        cmat = np.corrcoef(mat)
        cmat = np.nan_to_num(cmat)
        pca = PCA(n_components=1)
        mat_r = pca.fit(cmat).transform(cmat)
        eigs = np.array([t[0] for t in mat_r])
        eigs[ps] = 0
        if eig_r:  #flip the PC1 values according to other data such as histone markers
            eigs = 0 - eigs
    if virtual4C:
        predir = os.path.dirname(os.path.realpath(f))
        metaf = predir + "/petMeta.json"
        meta = json.loads(open(metaf).read())
        total = meta["Unique PETs"]
        #virtual4Csig = getVirtual4CSig(xy2, start, end, viewStart, viewEnd, res)
        virtual4Csig = getVirtual4CSig(xy2, start, end, viewStart, viewEnd)
        #virtual4Csig = virtual4Csig / total * 10**6
    if triu:
        mat = rotate(mat, angle=45, reshape=True)
        #take the uppper matrix and remove padding zeros
        to = int(mat.shape[0] / 2)
        mat = mat[to:, ]
        if norm == False:
            ns = list(range(mat.shape[0]))
            ns.reverse()
            for n in ns:
                if np.sum(mat[n, ]) > 0:
                    break
            mat = mat[:n, ]

            ns = list(range(mat.shape[1]))
            for na in ns:
                if np.sum(mat[:, na]) > 0:
                    break
            ns.reverse()
            for nb in ns:
                if np.sum(mat[:, nb]) > 0:
                    break
            mat = mat[:, na:nb + 1]
    #figure helights
    if triu:
        initSize = 3
        square = False
    else:
        initSize = 4
        square = False
    hights = initSize
    #heights ratio
    hr = []
    if gtf != "":
        genes = getGenes(gtf, chrom[0], start, end)
        """
        if len(genes) > 20:
            print(
                "More than 20 genes in the target region, only plot random 20."
            )
            ns = list(genes.keys())[:20]
            ng = {}
            for n in ns:
                ng[n] = genes[n]
            genes = ng
        """
        hights += len(genes) * 0.1
        hr.extend([0.1] * len(genes))
    if len(bws) > 0:
        hights += len(bws) * 0.5
        hr.extend([1] * len(bws))
    if oneD:
        hights += 0.5
        hr.append(1)
    if eig:
        hights += 0.5
        hr.append(1)
    if virtual4C:
        hights += 0.5
        hr.append(1)
    if loops is not None:
        hights += 0.5
        hr.append(1)
    if len(beds) > 0:
        hights += len(beds) * 0.2
        hr.extend([0.2] * len(beds))
    #heatmap and colorbar
    if triu:
        hr.extend([3, 0.1])
    else:
        hr.extend([6, 0.1])

    #prepare figure
    fig = pylab.figure(figsize=(width, hights))
    gs = mpl.gridspec.GridSpec(len(hr),
                               1,
                               height_ratios=hr,
                               top=0.95,
                               bottom=0.05,
                               left=0.1,
                               right=0.9,
                               wspace=0.05)
    pylab.suptitle(
        "%.2f kb,%s kb resolution, %s:%s-%s" %
        (float(end - start) / 1000.0, res / 1000.0, chrom[0], start, end),
        fontsize=8)
    axi = -1

    #plot gene
    if gtf != "":
        for n, g in genes.items():
            axi += 1
            ax = fig.add_subplot(gs[axi])
            plotGene(ax, n, g, start, end)

    #plot bigWig
    #prepare y-axis limitations
    bwvs = parseBwvs(bws, bwvs)
    #colors
    if bwcs == "":
        bwcs = range(len(bws))
    else:
        bwcs = list(map(int, bwcs.split(",")))
    for i, bw in enumerate(bws):
        axi += 1
        ax = fig.add_subplot(gs[axi])
        name = bw.split("/")[-1].split(".bw")[0]
        bw = pyBigWig.open(bw)
        ys = bw.values(chrom[0], start, end)
        ys = np.nan_to_num(ys)
        plotCoverage(ax,
                     ys,
                     colori=bwcs[i],
                     label=name,
                     vmin=bwvs[i][0],
                     vmax=bwvs[i][1])

    #plot 1D signal
    if oneD:
        axi += 1
        ax = fig.add_subplot(gs[axi])
        if oneDv != "":
            oneDv = list(map(float, oneDv.split(",")))
        else:
            oneDv = [None, None]
        plotCoverage(ax,
                     sig,
                     colori=3,
                     label="1D signal",
                     vmin=oneDv[0],
                     vmax=oneDv[1])

    #plot eigenvector
    if eig:
        axi += 1
        ax = fig.add_subplot(gs[axi])
        xs = np.arange(len(eigs))
        minxs = np.min(xs)
        maxns = np.max(xs)
        #ps = np.where( eigs!=0)[0]
        #eigs = eigs[ps]
        #xs = xs[ps]
        ps = np.where(eigs > 0)
        peigs = eigs[ps]
        pxs = xs[ps]
        ax.bar(pxs, peigs, color=colors[0], edgecolor=colors[0], alpha=0.5)
        ns = np.where(eigs < 0)
        neigs = eigs[ns]
        nxs = xs[ns]
        ax.bar(nxs, neigs, color=colors[1], edgecolor=colors[1], alpha=0.5)
        ax.plot(xs, [0] * len(xs), color="gray", alpha=0.8, linestyle="--")
        ax.tick_params(axis='both', which='major', labelsize=4)
        ax.set_xticklabels([])
        ax.set_xlim([minxs, maxns])
        ax.set_ylabel("eigenvector")

    #plot view point, virtual 4C plot
    if virtual4C:
        axi += 1
        ax = fig.add_subplot(gs[axi])
        if len(virtual4Csig) > 1000:
            virtual4Csig = getBinMean(virtual4Csig, 1000)
        #log2 is nesscessary
        virtual4Csig = np.log2(virtual4Csig + 1)
        xs = np.arange(len(virtual4Csig))
        ax.plot(xs, virtual4Csig, color=colors[0], label="virtual 4C signal")
        ax.fill_between(xs, 0, virtual4Csig, color=colors[0], alpha=0.8)
        ax.tick_params(axis='both', which='major', labelsize=4)
        ax.set_xticklabels([])
        ax.set_xlim([np.min(xs), np.max(xs)])
        ax.set_ylabel("log2(counts)", fontsize=6)
        ax.legend(fontsize=6, fancybox=False, frameon=False)

    #plot loops as arches
    nchrom = "-".join(chrom)
    if loops is not None and nchrom in loops and len(loops[nchrom]) > 0:
        axi += 1
        ax = fig.add_subplot(gs[axi])
        #plot the arch and annotate the PETs support the loop
        plotLoops(ax, loops, nchrom,start,end,xy2=xy2)
        
    #plot genomic features
    for i, bed in enumerate(beds):
        axi += 1
        name = bed.split("/")[-1].split(".bed")[0]
        ax = fig.add_subplot(gs[axi])
        rs = getBedRegion(bed, chrom[0], start, end)
        plotRegion(ax, rs, start, end, i, label=name)

    #plot the heatmap
    ax = fig.add_subplot(gs[-2])
    cax = fig.add_subplot(gs[-1])
    sns.set(font_scale=0.5)
    if corr:
        cmap = sns.diverging_palette(250, 15, s=75, l=40, n=11).as_hex()
        cmap[int(len(cmap) / 2)] = "#FFFFFF"
        cmap = ListedColormap(cmap)
        ax = sns.heatmap(mat,
                         xticklabels=False,
                         yticklabels=False,
                         square=square,
                         center=0,
                         linewidths=0.0,
                         ax=ax,
                         cmap=cmap,
                         cbar_ax=cax,
                         cbar_kws={
                             'label': ano,
                             'orientation': 'horizontal',
                             "shrink": 0.5,
                             "fraction": 0.2,
                             "anchor": (0.0, 1.0)
                         })

    else:
        if norm == False:
            #cmap = sns.cubehelix_palette(n_colors=11, as_cmap=True, light=1, hue=3)
            cmap = sns.light_palette("red", n_colors=9).as_hex()
            #cmap = plt.cm.Reds
            cmap[0] = "#FFFFFF"
            cmap = ListedColormap(cmap)
            center = None
            if vmin is None:
                vmin = 0
            vmax = vmax
        else:
            cmap = sns.color_palette("RdBu_r", 11).as_hex()
            cmap[int(len(cmap) / 2)] = "#FFFFFF"
            cmap = ListedColormap(cmap)
            center = 0
        ax = sns.heatmap(mat,
                         xticklabels=False,
                         yticklabels=False,
                         linewidths=0.0,
                         square=square,
                         center=center,
                         cmap=cmap,
                         vmin=vmin,
                         vmax=vmax,
                         ax=ax,
                         cbar_ax=cax,
                         cbar_kws={
                             'label': ano,
                             'orientation': 'horizontal',
                             "shrink": 0.3,
                             "fraction": 0.2,
                             "anchor": (0.0, 1.0)
                         })
    cax.tick_params(labelsize=4)

    #on the heatmap, draw the highlight region, such as TADs
    if domains != "":
        rs = getBedRegion(domains, chrom[0], start, end)
        if len(rs) > 0:
            pa = int(ax.get_xlim()[0])
            pb = int(ax.get_xlim()[1])
            if pa > pb:
                pa, pb = pb, pa
            ypa = int(ax.get_ylim()[0])
            ypb = int(ax.get_ylim()[1])
            if ypa > ypb:
                ypa, ypb = ypb, ypa
            for r in rs:
                npa = (r[0] - start) / (end - start)
                npb = (r[1] - start) / (end - start)
                if triu:
                    xa = npa * (pb - pa)
                    xb = npb * (pb - pa)
                    ya = npa * (ypb - ypa)
                    yb = npb * (ypb - ypa)
                    ax.plot([xa, (xa + xb) / 2, xb], [ypa, (yb - ya), ypa],
                            color=colors[1],
                            linewidth=1,
                            linestyle="--")
                else:
                    ax.axvline(x=npa * (pb - pa),
                               ymin=1 - npa,
                               ymax=1 - npb,
                               color=colors[1],
                               linewidth=1,
                               linestyle="--")
                    ax.axvline(x=npb * (pb - pa),
                               ymin=1 - npa,
                               ymax=1 - npb,
                               color=colors[1],
                               linewidth=1,
                               linestyle="--")
                    ax.axhline(y=npa * (pb - pa),
                               xmin=npa,
                               xmax=npb,
                               color=colors[1],
                               linewidth=1,
                               linestyle="--")
                    ax.axhline(y=npb * (pb - pa),
                               xmin=npa,
                               xmax=npb,
                               color=colors[1],
                               linewidth=1,
                               linestyle="--")
    #draw the box
    ax.axvline(x=ax.get_xlim()[0], color="k", linewidth=2)
    ax.axvline(x=ax.get_xlim()[1], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[0], color="k", linewidth=2)
    ax.axhline(y=ax.get_ylim()[1], color="k", linewidth=2)
    if not triu:
        pylab.tight_layout()
    pylab.savefig(fo + "_matrix.pdf")


def plotPETsArches(
        f,
        fo,
        start=0,
        end=-1,
        cut=0,
        mcut=-1,
        oneD=False,
        oneDv="",
        bws=[],
        bwvs="",
        bwcs="",
        beds=[],
        loops=None,
        gtf="",
        aw=1,
        ac=1,
        aa=1,
        width=4,
):
    """
    Plot the interacting PETs as arches, showing the raw data. 
    """
    chrom, xy = parseIxy(f, cut=cut, mcut=mcut)
    xy2 = XY(xy[:, 0], xy[:, 1])  #XY object

    if start == 0:
        start = np.min(xy)
    if end == -1:
        end = np.max(xy)
    if oneD:
        predir = os.path.dirname(os.path.realpath(f))
        metaf = predir + "/petMeta.json"
        meta = json.loads(open(metaf).read())
        total = meta["Unique PETs"] * 2
        sig = get1DSig(xy2, start, end)
        sig = sig / total * 10**6
    hights = 0
    #heights ratio
    hr = []
    if gtf != "":
        genes = getGenes(gtf, chrom[0], start, end)
        """
        if len(genes) > 20:
            print(
                "More than 20 genes in the target region, only plot random 20."
            )
            ns = list(genes.keys())[:20]
            ng = {}
            for n in ns:
                ng[n] = genes[n]
            genes = ng
        """
        hights += len(genes) * 0.1
        hr.extend([0.1] * len(genes))
    if len(bws) > 0:
        hights += len(bws) * 0.5
        hr.extend([1] * len(bws))
    if oneD:
        hights += 0.5
        hr.append(1)
    if loops is not None:
        hights += 0.5
        hr.append(1)
    if len(beds) > 0:
        hights += len(beds) * 0.2
        hr.extend([0.2] * len(beds))
    #arches
    hights += 2
    hr.append(2.5)

    #prepare figure
    fig = pylab.figure(figsize=(width, hights))
    gs = mpl.gridspec.GridSpec(len(hr),
                               1,
                               height_ratios=hr,
                               top=0.9,
                               bottom=0.05,
                               left=0.1,
                               right=0.9,
                               wspace=0.05)
    pylab.suptitle("%.2f kb,%s:%s-%s" %
                   (float(end - start) / 1000.0, chrom[0], start, end),
                   fontsize=8)
    axi = -1

    #plot gene
    if gtf != "":
        for n, g in genes.items():
            axi += 1
            ax = fig.add_subplot(gs[axi])
            plotGene(ax, n, g, start, end)

    #plot bigWig
    #yaxis limitaitons
    bwvs = parseBwvs(bws, bwvs)
    #colors
    if bwcs == "":
        bwcs = range(len(bws))
    else:
        bwcs = list(map(int, bwcs.split(",")))
    for i, bw in enumerate(bws):
        axi += 1
        ax = fig.add_subplot(gs[axi])
        name = bw.split("/")[-1].split(".bw")[0]
        bw = pyBigWig.open(bw)
        ys = bw.values(chrom[0], start, end)
        ys = np.nan_to_num(ys)
        plotCoverage(ax,
                     ys,
                     colori=bwcs[i],
                     label=name,
                     vmin=bwvs[i][0],
                     vmax=bwvs[i][1])

    #plot 1D signal
    if oneD:
        axi += 1
        ax = fig.add_subplot(gs[axi])
        if oneDv != "":
            oneDv = list(map(float, oneDv.split(",")))
        else:
            oneDv = [None, None]
        plotCoverage(ax,
                     sig,
                     colori=3,
                     label="1D signal",
                     vmin=oneDv[0],
                     vmax=oneDv[1])

    #plot loops as arches
    nchrom = "-".join(chrom)
    if loops is not None and nchrom in loops and len(loops[nchrom]) > 0:
        axi += 1
        ax = fig.add_subplot(gs[axi])
        #plot the arch and annotate the PETs support the loop
        #get the minal PETs number as linewidth 1,others are fold
        plotLoops(ax, loops, nchrom,start,end,xy2=xy2)
       
    #plot genomic features
    for i, bed in enumerate(beds):
        axi += 1
        name = bed.split("/")[-1].split(".bed")[0]
        ax = fig.add_subplot(gs[axi])
        rs = getBedRegion(bed, chrom[0], start, end)
        plotRegion(ax, rs, start, end, i, label=name)

    #plot PETs as arches
    axi += 1
    ax = fig.add_subplot(gs[axi])
    ps = xy2.queryPeakBoth(start, end)
    pa = 0
    pb = 1.0
    if len(ps) > 0:
        ymax = 0
        for p in ps:
            ca = xy[p, 0]
            cb = xy[p, 1]
            cc = (ca + cb) / 2
            npa = float(ca - start) / (end - start) * (pb - pa)
            npb = float(cb - start) / (end - start) * (pb - pa)
            npc = float(cc - start) / (end - start) * (pb - pa)
            a = npb - npa  #a is x axis size for eclipse
            b = a / 2 * 0.6  #b is y axis size for eclipse
            if b > ymax:
                ymax = b
            ax.add_patch(
                Arc(
                    (npc, 0),
                    a,
                    b,
                    theta1=180,
                    theta2=360,
                    edgecolor=colors[ac],
                    lw=aw,
                    alpha=aa,
                ))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim([0, 1])
        ax.set_ylim([0, -ymax * 0.55])
        ax.invert_yaxis()
    pylab.savefig(fo + "_arches.pdf")


def plotPETsScatter(
        f,
        fo,
        start=0,
        end=-1,
        cut=0,
        mcut=-1,
        oneD=False,
        oneDv="",
        bws=[],
        bwvs="",
        bwcs="",
        beds=[],
        loops=None,
        gtf="",
        ss = 1,
        sc = 0,
        sa = 0.5,
        triu=False,
        width=8,
):
    """
    Plot the interacting PETs as scatter, showing the raw data. 
    """
    chrom, xy = parseIxy(f, cut=cut, mcut=mcut)
    xy2 = XY(xy[:, 0], xy[:, 1])  #XY object

    if start == 0:
        start = np.min(xy)
    if end == -1:
        end = np.max(xy)
    if oneD:
        predir = os.path.dirname(os.path.realpath(f))
        metaf = predir + "/petMeta.json"
        meta = json.loads(open(metaf).read())
        total = meta["Unique PETs"] * 2
        sig = get1DSig(xy2, start, end)
        sig = sig / total * 10**6
    hights = 0
    #heights ratio
    hr = []
    if gtf != "":
        genes = getGenes(gtf, chrom[0], start, end)
        hights += len(genes) * 0.1
        hr.extend([0.1] * len(genes))
    if len(bws) > 0:
        hights += len(bws) * 0.5
        hr.extend([1] * len(bws))
    if oneD:
        hights += 0.5
        hr.append(1)
    if loops is not None:
        hights += 0.5
        hr.append(1)
    if len(beds) > 0:
        hights += len(beds) * 0.2
        hr.extend([0.2] * len(beds))
    #scatter plot
    hights += 2
    hr.append(2.5)
    if triu:
        hights += 3
        square = False
    else:
        hights += 4
        square = False
  
    #prepare figure
    fig = pylab.figure(figsize=(width, hights))
    gs = mpl.gridspec.GridSpec(len(hr),
                               1,
                               height_ratios=hr,
                               top=0.9,
                               bottom=0.05,
                               left=0.1,
                               right=0.9,
                               wspace=0.05)
    pylab.suptitle("%.2f kb,%s:%s-%s" %
                   (float(end - start) / 1000.0, chrom[0], start, end),
                   fontsize=8)
    axi = -1

    #plot gene
    if gtf != "":
        for n, g in genes.items():
            axi += 1
            ax = fig.add_subplot(gs[axi])
            plotGene(ax, n, g, start, end)

    #plot bigWig
    #yaxis limitaitons
    bwvs = parseBwvs(bws, bwvs)
    #colors
    if bwcs == "":
        bwcs = range(len(bws))
    else:
        bwcs = list(map(int, bwcs.split(",")))
    for i, bw in enumerate(bws):
        axi += 1
        ax = fig.add_subplot(gs[axi])
        name = bw.split("/")[-1].split(".bw")[0]
        bw = pyBigWig.open(bw)
        ys = bw.values(chrom[0], start, end)
        ys = np.nan_to_num(ys)
        plotCoverage(ax,
                     ys,
                     colori=bwcs[i],
                     label=name,
                     vmin=bwvs[i][0],
                     vmax=bwvs[i][1])

    #plot 1D signal
    if oneD:
        axi += 1
        ax = fig.add_subplot(gs[axi])
        if oneDv != "":
            oneDv = list(map(float, oneDv.split(",")))
        else:
            oneDv = [None, None]
        plotCoverage(ax,
                     sig,
                     colori=3,
                     label="1D signal",
                     vmin=oneDv[0],
                     vmax=oneDv[1])

    #plot loops as arches
    nchrom = "-".join(chrom)
    if loops is not None and nchrom in loops and len(loops[nchrom]) > 0:
        axi += 1
        ax = fig.add_subplot(gs[axi])
        #plot the arch and annotate the PETs support the loop
        #get the minal PETs number as linewidth 1,others are fold
        plotLoops(ax, loops, nchrom,start,end,xy2=xy2)
       
    #plot genomic features
    for i, bed in enumerate(beds):
        axi += 1
        name = bed.split("/")[-1].split(".bed")[0]
        ax = fig.add_subplot(gs[axi])
        rs = getBedRegion(bed, chrom[0], start, end)
        plotRegion(ax, rs, start, end, i, label=name)

    #plot PETs as dos
    axi += 1
    ax = fig.add_subplot(gs[axi])
    ps = xy2.queryPeakBoth(start, end)
    mat = xy[list(ps)] 
    mat = mat - start
    if triu:
        #caculating the rotate coordinates
        x = mat[:,0]*np.cos( -np.pi/4 ) - mat[:,1]*np.sin( -np.pi/4 )
        y = mat[:,1]*np.cos( -np.pi/4 ) + mat[:,0]*np.sin( -np.pi/4 )
        ax.set_xlim([0, np.max(x)])
        ax.set_ylim([np.min(y),np.max(y)])
        ax.scatter( x,y, s =ss, color=colors[sc], alpha=sa)
    else:
        ax.scatter( mat[:,0], mat[:,1], s =ss, color=colors[sc], alpha=sa)
        ax.scatter( mat[:,1], mat[:,0], s =ss, color=colors[sc], alpha=sa)
        ax.set_xlim([0, end-start])
        ax.set_ylim([0, end-start])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    pylab.savefig(fo + "_scatter.pdf")


def plotProfiles(
        fo,
        chrom="",
        start=0,
        end=-1,
        bws=[],
        bwvs="",
        bwcs="",
        beds=[],
        loops=None,
        gtf="",
        width=8,
):
    """
    Plot profiles. 
    """
    #heights ratio
    hights = 0
    hr = []
    if gtf != "":
        genes = getGenes(gtf, chrom, start, end)
        """
        if len(genes) > 20:
            print(
                "More than 20 genes in the target region, only plot random 20."
            )
            ns = list(genes.keys())[:20]
            ng = {}
            for n in ns:
                ng[n] = genes[n]
            genes = ng
        """
        hights += len(genes) * 0.12
        hr.extend([0.12] * len(genes))
    if len(bws) > 0:
        hights += len(bws) * 0.3
        hr.extend([0.8] * len(bws))
    if loops is not None:
        hights += 0.5
        hr.append(1)
    if len(beds) > 0:
        hights += len(beds) * 0.2
        hr.extend([0.2] * len(beds))

    #prepare figure
    fig = pylab.figure(figsize=(width, hights))
    gs = mpl.gridspec.GridSpec(
        len(hr),
        1,
        height_ratios=hr,
        top=0.9,
        bottom=0.05,
        left=0.1,
        right=0.9,
        wspace=0.0,
        hspace=0.05,
    )
    pylab.suptitle("%.2f kb,%s:%s-%s" %
                   (float(end - start) / 1000.0, chrom, start, end),
                   fontsize=8)
    axi = -1
    #plot gene
    if gtf != "":
        for n, g in genes.items():
            axi += 1
            ax = fig.add_subplot(gs[axi])
            plotGene(ax, n, g, start, end)

    #plot bigWig
    #yaxis limitaitons
    bwvs = parseBwvs(bws, bwvs)
    #colors
    if bwcs == "":
        bwcs = range(len(bws))
    else:
        bwcs = list(map(int, bwcs.split(",")))
    for i, bw in enumerate(bws):
        axi += 1
        ax = fig.add_subplot(gs[axi])
        name = bw.split("/")[-1].split(".bw")[0]
        bw = pyBigWig.open(bw)
        ys = bw.values(chrom, start, end)
        ys = np.nan_to_num(ys)
        ax = plotCoverage(ax,
                          ys,
                          colori=bwcs[i],
                          label=name,
                          vmin=bwvs[i][0],
                          vmax=bwvs[i][1])
        if i == 0:
            sns.despine(ax=ax, bottom=False, right=False, left=False, top=False)
        elif i == len(bws) - 1:
            sns.despine(ax=ax, bottom=False, right=False, left=False, top=True)
        else:
            #sns.despine(ax=ax, bottom=True, right=False, left=False, top=True)
            sns.despine(ax=ax, bottom=False, right=False, left=False, top=True)

    #plot loops as arches
    nchrom = chrom + "-" + chrom
    if loops is not None and nchrom in loops and len(loops[nchrom]) > 0:
        axi += 1
        ax = fig.add_subplot(gs[axi])
        #plot the arch for loops, all same width
        plotLoops(ax, loops, nchrom,start,end)

    #plot genomic features
    for i, bed in enumerate(beds):
        axi += 1
        name = bed.split("/")[-1].split(".bed")[0]
        ax = fig.add_subplot(gs[axi])
        rs = getBedRegion(bed, chrom, start, end)
        plotRegion(ax, rs, start, end, i, label=name)
    pylab.savefig(fo + "_profiles.pdf")
