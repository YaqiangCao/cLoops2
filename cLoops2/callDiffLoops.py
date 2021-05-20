#!/usr/bin/env python
#--coding:utf-8 --
"""
callDiffLoops.py
Calling differentially enriched loops between different cells/conditions or similar situation. 
Well tested for Trac-looping data.  
2020-03-10: add MA and voccano plot.
2020-05-04: MANorm2 normalization added.
2020-05-05: refined a lot with the base of MANorm
2020-05-06: updated auto-estimation of MA Mcut and Acut
2020-05-10: fitting values changed to counts , if using interaction per kb, there will be outlieers in MA plot. The cutoffs of M and A is determine in background, while there is systematic difference between bg and fg, a ftting again is needed.
2020-05-20: try to add the shift of background to foreground, too strigenet
2020-06-24: try to add the estimation of anchors
2020-07-30: try to add the estimation of anchors, and seperated loops. acut and mcut estimated from background already very strong. win size may not affect, so for efficient consideration set to 1. Alsoadd p-values to bg estimation.
2021-04-12: cutomize parameters for acut and mcut for MA plot added
2021-05-10: add heatmap vmin/vmax, cmap; not p-value cutoff added to include more loops.
"""

__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os
import sys
import json
from glob import glob
from copy import deepcopy
from datetime import datetime

#3rd
import joblib
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import poisson
from sklearn import linear_model
from joblib import Parallel, delayed
from sklearn.mixture import GaussianMixture as GMM

#cLoops
from cLoops2.io import parseTxt2Loops, ixy2pet, dloops2txt, dloops2juiceTxt, loops2washuTxt, dloops2NewWashuTxt
from cLoops2.ds import Loop, XY, DiffLoop
from cLoops2.agg import getALoops
from cLoops2.geo import checkLoopOverlap
from cLoops2.est import estSfMANorm
from cLoops2.settings import *
from cLoops2.callCisLoops import getLoopNearbyPETs


def mergeLoops(aloops, samplea, bloops, sampleb):
    """
    Get the union set of loops for candidate test samples. 
    """
    mloops = {}
    allkeys = set(aloops.keys())
    allkeys.update(bloops.keys())
    for key in allkeys:
        mloops[key] = []

    for key in mloops.keys():
        if key not in aloops:
            for loop in bloops[key]:
                loop.id = loop.id + "|" + sampleb
                mloops[key].append(loop)
        elif key not in bloops:
            for loop in aloops[key]:
                loop.id = loop.id + "|" + samplea
                mloops[key].append(loop)
        else:
            nloops = []
            #scan a first
            for aloop in aloops[key]:
                flag = True  # no overlapped loops, keep it
                p = -1
                for j, bloop in enumerate(bloops[key]):
                    if checkLoopOverlap(aloop, bloop):
                        flag = False
                        p = j
                #no overlapped loops find, just add it
                if flag == True:
                    aloop.id = aloop.id + "|" + samplea
                    nloops.append(aloop)
                else:  #merge overlapped loops, skip the bloop
                    bloop = deepcopy(bloops[key][p])
                    del bloops[key][p]
                    nloop = Loop()
                    nloop.id = "|".join([aloop.id, samplea, bloop.id, sampleb])
                    nloop.chromX = aloop.chromX
                    nloop.chromY = aloop.chromY
                    nloop.x_start = min(aloop.x_start, bloop.x_start)
                    nloop.x_end = max(aloop.x_end, bloop.x_end)
                    nloop.x_center = (nloop.x_start + nloop.x_end) / 2
                    nloop.y_start = min(aloop.y_start, bloop.y_start)
                    nloop.y_end = max(aloop.y_end, bloop.y_end)
                    nloop.y_center = (nloop.y_start + nloop.y_end) / 2
                    if set(key.split("-")) == 1:
                        nloop.distance = nloop.y_center - nloop.x_center
                    nloops.append(nloop)
            #scan b to add to merged
            for bloop in bloops[key]:
                bloop.id = bloop.id + "|" + sampleb
                nloops.append(bloop)
            mloops[key] = nloops
    return mloops


def quantDloops(key, loops, fixya, fixyb, cut=0, mcut=-1, win=1, pseudo=1.0):
    """
    Estimate the differential of significant loops.
    """
    if len(loops) == 0:
        return None

    #query structure
    axy = ixy2pet(fixya, cut=cut, mcut=mcut)
    bxy = ixy2pet(fixyb, cut=cut, mcut=mcut)

    #all information of nearby counts
    ts = []
    cs = []

    print("Quantify %s loops for %s." % (len(loops), key))
    #for loop in tqdm(loops):
    dloops = []
    for loop in tqdm(loops):
        ara, arb, arab = axy.queryLoop(loop.x_start, loop.x_end, loop.y_start,
                                       loop.y_end)
        bra, brb, brab = bxy.queryLoop(loop.x_start, loop.x_end, loop.y_start,
                                       loop.y_end)
        arabs, anbps = getLoopNearbyPETs(loop, axy, win)
        brabs, bnbps = getLoopNearbyPETs(loop, bxy, win)
        #get all the information
        dloop = DiffLoop()
        dloop.id = loop.id
        dloop.chromX = loop.chromX
        dloop.chromY = loop.chromY
        dloop.x_start = loop.x_start
        dloop.x_end = loop.x_end
        dloop.x_center = loop.x_center
        dloop.y_start = loop.y_start
        dloop.y_end = loop.y_end
        dloop.y_center = loop.y_center
        dloop.distance = loop.y_center - loop.x_center
        dloop.raw_trt_ra = len(ara)
        dloop.raw_trt_rb = len(arb)
        dloop.raw_con_ra = len(bra)
        dloop.raw_con_rb = len(brb)
        dloop.raw_trt_rab = len(arab)
        dloop.raw_con_rab = len(brab)
        dloop.size = dloop.x_end - dloop.x_start + dloop.y_end - dloop.y_start
        dloop.raw_trt_mrab = np.mean(arabs)
        dloop.raw_con_mrab = np.mean(brabs)
        ts.extend(arabs)
        cs.extend(brabs)
        dloop.trt_es = dloop.raw_trt_rab / max(dloop.raw_trt_mrab, pseudo)
        dloop.con_es = dloop.raw_con_rab / max(dloop.raw_con_mrab, pseudo)
        dloops.append(dloop)
    return key, ts, cs, dloops


def plotBgMANorm(sf, cs, ts, m, a, m2, a2, fout, mcut, acut, fdr):
    """
    Plot the MANorm2 result before and after for background data, and estimate cutoffs
    @param sf: [], list of sacling factors, 
    @param cs: [], list of control background data,log2 transformed
    @param ts: [], list of target background data,log2 transformed
    @param m: [], np.log2(cs)-np.log2(ts)
    @param a: [], (np.log2(cs)+np.log2(ts))/2
    """
    fig, axs = pylab.subplots(1,
                              3,
                              figsize=(8.5, 2.2),
                              sharex=False,
                              sharey=False)
    axs = axs.reshape(-1)
    #plot 1 raw data
    axs[0].scatter(ts, cs, s=1, c="gray")
    axs[0].set_xlabel("target log2(PETs)")
    axs[0].set_ylabel("control log2(PETs)")
    axs[0].set_title("raw bg data")
    #plot the fit
    #show the formula
    if sf[1] > 0:
        formula = "y=%.3fx+%.3f" % (sf[0], sf[1])
    else:
        formula = "y=%.3fx%.3f" % (sf[0], sf[1])
    axs[0].plot(ts, ts * sf[0] + sf[1], label=formula)
    axs[0].legend()
    #plot  2, raw data MA
    axs[1].scatter(a, m, s=1, c="gray")
    axs[1].set_title("before normalization\nmean(M):%.3f;M~A PCC:%.3f" %
                     (np.mean(m), np.corrcoef(m, a)[0][1]),
                     fontsize=8)
    axs[1].set_xlabel("A")
    axs[1].set_ylabel("M")
    #plot 3 transformed data
    axs[2].scatter(a2, m2, s=1, c="gray")
    axs[2].set_title("after normalization\nmean(M):%.3f;M~A PCC :%.3f" %
                     (np.mean(m2), np.corrcoef(m2, a2)[0][1]),
                     fontsize=8)

    axs[2].set_xlabel("A")
    axs[2].set_ylabel("M")
    upm = np.where(m2 >= mcut)[0]
    upa = np.where(a2 >= acut)[0]
    up = list(set(upm).intersection(upa))
    downm = np.where(m2 <= -mcut)[0]
    downa = np.where(a2 >= acut)[0]
    down = list(set(downm).intersection(downa))
    #up dots
    axs[2].scatter(a2[up],
                   m2[up],
                   color=colors[2],
                   s=1,
                   label="up %.3f%s" % (float(len(up)) / len(m2) * 100, "%"))
    #down dots
    axs[2].scatter(a2[down],
                   m2[down],
                   color=colors[3],
                   s=1,
                   label="down %.3f%s" %
                   (float(len(down)) / len(m2) * 100, "%"))
    axs[2].legend()
    pylab.suptitle("bg FDR<=%.3f mcut=%.3f acut=%.3f" % (fdr, mcut, acut),
                   fontsize=8)
    pylab.subplots_adjust(top=0.7, wspace=0.3, hspace=0.5)
    pylab.savefig("%s_background_MANormFit.pdf" % fout)


def getBgNorm(cs,
              ts,
              fout,
              mcut=0,
              acut=0,
              step=0.1,
              fdrcut=0.05,
              pseudocut=1.0
    ):
    """ 
    Do the MANorm with background data and estimate the cutoffs.
    @param cs: [], list of background counts of control data
    @param ts: [], list of background counts of target data
    @param cfgs: [], list of foreground counts of control data
    @param tfgs: [], list of foreground counts of target data
    @param fout: output prefix for plot
    @param mcut: float,M cutoff for MA plot
    @param acut: float,A cutoff for MA plot
    @param step: float, step for increasing mcut and acut while searching
    @param fdrcut: float, main cutoff
    @param pseudocut: float, used to avoid 0 
    """
    #remove zeros for log transformation
    cs = np.array(cs)
    ts = np.array(ts)
    cs = np.log2(cs + pseudocut)
    ts = np.log2(ts + pseudocut)

    #linear shiftting
    sf = estSfMANorm(cs, ts)
    m = cs - ts
    a = (ts + cs) / 2
    #transform the data
    ts2 = [sf[0] * t + sf[1] for t in ts]

    m2 = cs - ts2
    a2 = (cs + ts2) / 2
    m2abs = np.abs(m2)
    #get the acut
    fdr = 1
    i = 0
    while fdr > fdrcut:
        if i % 2 == 0:
            acut = acut + step
        else:
            mcut = mcut + step
        i += 1
        de = np.where(m2abs >= mcut)[0]
        de = len(np.where(a2[de] >= acut)[0])
        fdr = float(de) / len(m2abs)
    #plot the estimated result
    plotBgMANorm(sf, cs, ts, m, a, m2, a2, fout, mcut, acut, fdr)
    return sf, acut, mcut


def estLoopDiffSig(key, sf, ta, tb, dloops, pseudo=1.0):
    """
    Estimation of differential significance.
    @param key: str,chrom-chrom
    @param sf: [float,float],scaling factor,
    @param ta: float or int, total PETs for sample A
    @param tb: float or int, total PETs for sample B
    @param dloops: list of cLoops2.ds.DiffLoop object
    @param pseudo: float/int, pseudo counts, used to avoid /0 or log(0)
    """
    #cs for control ipk, ts for target ipk, ts2 for scaled
    cs, ts, ts2 = [], [], []
    print("Estimate difference significance %s loops for %s." %
          (len(dloops), key))
    for loop in tqdm(dloops):
        loop.trt_density = max(loop.raw_trt_rab,
                               pseudo) / loop.size / ta * 10**9
        loop.con_density = max(loop.raw_con_rab,
                               pseudo) / loop.size / tb * 10**9
        #the MA fitting is based on log2 transformation
        loop.scaled_trt_rab = 2**(sf[0] * np.log2(loop.raw_trt_rab) + sf[1])
        loop.scaled_trt_mrab = 2**(sf[0] * np.log2(loop.raw_trt_mrab) + sf[1])
        loop.scaled_trt_ra = 2**(sf[0] * np.log2(loop.raw_trt_ra) + sf[1])
        loop.scaled_trt_rb = 2**(sf[0] * np.log2(loop.raw_trt_rb) + sf[1])
        #target sample has potential to be significant
        if loop.scaled_trt_rab > loop.raw_con_rab:
            fg = loop.scaled_trt_rab
            bg = max(loop.raw_con_rab, loop.raw_con_mrab)
        #control sample has potential to be significant
        else:
            fg = loop.raw_con_rab
            bg = max(loop.scaled_trt_rab, loop.scaled_trt_mrab)
        cs.append(max(loop.raw_con_rab, pseudo))
        ts.append(max(loop.raw_trt_rab, pseudo))
        ts2.append(max(loop.scaled_trt_rab, pseudo))
        pop = poisson.sf(fg - 1, bg)
        pop = max([pop, 1e-300])
        loop.poisson_p_value = pop
        #loop.raw_fc = np.log2(max(loop.raw_trt_rab, pseudo) / ta) - np.log2( max(loop.raw_con_rab, pseudo) / tb)
        #loop.scaled_fc = np.log2(max(loop.scaled_trt_rab, pseudo)) - np.log2( max(loop.raw_con_rab, pseudo) )
        loop.raw_fc = np.log2(loop.raw_trt_rab / ta) - np.log2(
            loop.raw_con_rab / tb)
        loop.scaled_fc = np.log2(loop.scaled_trt_rab) - np.log2(
            loop.raw_con_rab)
    return dloops, cs, ts, ts2


def markDiffSig(loops, acut, mcut, pcut=1e-2, pseudo=1,igp=False,noPCorr=False):
    """
    Carry out Bonferroni correction for p-values first then mark the significance of loops
    """
    for loop in loops:
        if noPCorr == False:
            loop.poisson_p_value = min(1, loop.poisson_p_value * len(loops))
        if igp == False:
            if loop.poisson_p_value <= pcut and abs(loop.scaled_fc) >= mcut:
                c = np.log2(max(loop.raw_con_rab, pseudo))
                t = np.log2(max(loop.scaled_trt_rab, pseudo))
                #loop.significant = 1
                a = (c + t) / 2
                if a >= acut:
                    loop.significant = 1
                else:
                    loop.significant = 0
            else:
                loop.significant = 0
        else:
            if abs(loop.scaled_fc) >= mcut:
                c = np.log2(max(loop.raw_con_rab, pseudo))
                t = np.log2(max(loop.scaled_trt_rab, pseudo))
                #loop.significant = 1
                a = (c + t) / 2
                if a >= acut:
                    loop.significant = 1
                else:
                    loop.significant = 0
            else:
                loop.significant = 0
    return loops


def plotDiffLoopsMA(sigIndex, cs, ts, ts2, tname, cname, mcut, acut, fout):
    """
    Plot the MA plot for differential enriched loops.
    """
    cs = np.log2(cs)
    ts = np.log2(ts)
    ts2 = np.log2(ts2)
    m = ts - cs
    a = (ts + cs) / 2
    m2 = ts2 - cs
    a2 = (ts2 + cs) / 2
    #start plot
    fig, axs = pylab.subplots(1, 2, figsize=(6.4, 2.75 * 0.8))
    #raw data
    axs[0].scatter(a, m, s=1, c="gray")
    axs[0].set_title("before normalization\nmean(M):%.3f;M~A PCC:%.3f" %
                     (np.mean(m), np.corrcoef(m, a)[0][1]),
                     fontsize=8)
    axs[0].set_xlabel("A, 1/2( log2(%s)+log2(%s) )" % (tname, cname),
                      fontsize=6)
    axs[0].set_ylabel("M, log2(%s) - log2(%s)" % (tname, cname), fontsize=6)
    #scaled data
    axs[1].scatter(a2, m2, s=1, c="gray")
    axs[1].set_title("after normalization\nmean(M):%.3f;M~A PCC :%.3f" %
                     (np.mean(m2), np.corrcoef(m2, a2)[0][1]),
                     fontsize=8)
    #diffrentially enriched loops
    up = np.where(m2 >= mcut)[0]
    up = list(set(up).intersection(set(sigIndex)))
    up = list(set(np.where(a2 >= acut)[0]).intersection(up))
    axs[1].scatter(a2[up], m2[up], color=colors[2], s=1, alpha=1)
    down = np.where(m2 <= -mcut)[0]
    down = list(set(down).intersection(set(sigIndex)))
    down = list(set(np.where(a2 >= acut)[0]).intersection(down))
    #cutoff lines
    axs[1].scatter(a2[down], m2[down], color=colors[3], s=2, alpha=1)
    axs[1].axhline(y=0,
                   linewidth=1,
                   linestyle="--",
                   color=colors[0],
                   alpha=0.5)
    axs[1].axhline(y=mcut, linewidth=1, linestyle="--", color=colors[1])
    axs[1].axhline(y=-mcut, linewidth=1, linestyle="--", color=colors[1])
    axs[1].axvline(x=acut,
                   linewidth=1,
                   linestyle="--",
                   color=colors[4],
                   alpha=0.5)
    mm = np.max(m2) * 0.8
    ma = np.max(a2) * 0.7
    axs[1].text(ma, mm, "%s loops" % len(up), color=colors[2])
    axs[1].text(ma, -mm, "%s loops" % len(down), color=colors[3])
    axs[1].set_xlabel("A, 1/2( log2(%s)+log2(%s) )" % (tname, cname),
                      fontsize=6)
    axs[1].set_ylabel("M, log2(%s) - log2(%s)" % (tname, cname), fontsize=6)
    fig.tight_layout()
    pylab.savefig("%s_diffLoopsMA.pdf" % (fout))


def plotDiffLoopsVolcano(f, output, tname, cname, fccut=1, pcut=1e-2):
    """
    Plot the MA plot for differential enriched loops.
    """
    mat = pd.read_csv(f, index_col=0, sep="\t")
    fig, ax = pylab.subplots()
    fc = mat["scaledFc"]
    ps = mat["poissonPvalue"]
    ps = -np.log10(ps)
    s = mat["significant"]
    s = s[s > 0].index
    #all dots
    ax.scatter(fc, ps, color="gray", s=1, alpha=0.5)
    up = fc[fc > 0].index.intersection(s)
    down = fc[fc < 0].index.intersection(s)
    #up dots
    ax.scatter(fc[up], ps[up], color=colors[2], s=1)
    #down dots
    ax.scatter(fc[down], ps[down], color=colors[3], s=1)
    ax.axhline(y=-np.log10(pcut), linewidth=1, linestyle="--", color=colors[0])

    ax.text(3, 90, "%s loops" % len(up), color=colors[2])
    ax.text(-3, 90, "%s loops" % len(down), color=colors[3])

    ax.set_xlabel("log2( %s/%s )" % (tname, cname))
    ax.set_ylabel("-log10(p-value)")
    ax.set_ylim([-1, 100])
    pylab.savefig("%s_diffLoopsVolcano.pdf" % (output))


def getDiffAggLoops(predir, loops, cpu=1, norm=True):
    """
    Get the mean matrix and enrichment score.
    """
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(loops.keys())))
    ds = Parallel(n_jobs=cpu, backend="multiprocessing")(delayed(getALoops)(
        key,
        meta["data"]["cis"][key]["ixy"],
        loops[key],
    ) for key in keys)
    mat = np.concatenate([d[0] for d in ds if d[0] is not None], axis=0)
    es = []
    for i in range(mat.shape[0]):
        p = int(mat[i].shape[0] / 2)
        if mat[i].mean() > 0:
            nmat = deepcopy(mat[i])
            nmat[p, p] = 0
            if nmat.mean() == 0.0:
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
    return mat, es


def plotDiffAggLoops(dloops, output, tl, cl, td, cd, cpu=1, norm=True,vmin=None,vmax=None,cmap="summer"):
    """
    Plot the aggregated unqiue and overlapped loops.
    """
    #process meta information
    na = td.split("/")[-1]  #name of sample directory
    nb = cd.split("/")[-1]

    #seperated loops as overlapped, trt specific, con specific
    overlappedLoops = {}
    trtLoops = {}
    conLoops = {}
    #counts of called un-specific and specific loops
    cover, ctrt, ccon = 0, 0, 0
    for loop in dloops:
        key = loop.chromX + "-" + loop.chromY
        if loop.significant < 1:
            if key not in overlappedLoops:
                overlappedLoops[key] = []
            overlappedLoops[key].append(loop)
            cover += 1
        else:
            if loop.scaled_fc > 0:
                if key not in trtLoops:
                    trtLoops[key] = []
                trtLoops[key].append(loop)
                ctrt += 1
            else:
                if key not in conLoops:
                    conLoops[key] = []
                ccon += 1
                conLoops[key].append(loop)
    
    #cmap for heatmaps
    if cmap is None:
        cmap = cmaps["summer"]
    else:
        cmap = cmaps[cmap]

    #show enrichment score and aggregation plot
    fig, axs = pylab.subplots(2, 3, figsize=(10, 6))
    
    ax = axs[0][0]
    if cover > 0:
        trtOverMat, trtOverES = getDiffAggLoops(td, overlappedLoops, cpu)
        sns.heatmap(trtOverMat,
                    xticklabels=False,
                    yticklabels=False,
                    square=True,
                    ax=ax,
                    cmap=cmap,
                    linewidths=0.05,
                    linecolor="gray",
                    linestyle="--",
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={"shrink": 0.5})
        ax.set_ylabel(na, fontsize=10)
        ax.set_title("%s un-specific loops\nES:%.3f" %
                     (cover, np.mean(trtOverES)),
                     fontsize=8)
    else:
        ax.set_title("No common loops")

    ax = axs[0][1]
    if ctrt > 0:
        trtTrtMat, trtTrtES = getDiffAggLoops(td, trtLoops, cpu)
        sns.heatmap(trtTrtMat,
                    xticklabels=False,
                    yticklabels=False,
                    square=True,
                    ax=ax,
                    cmap=cmap,
                    linewidths=0.05,
                    linecolor="gray",
                    linestyle="--",
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={"shrink": 0.5})
        ax.set_title("%s specific loops\nES:%.3f" %
                     (ctrt, np.mean(trtTrtES)),
                     fontsize=8)
    else:
        ax.set_title("No %s unique loops" % na)

    ax = axs[0][2]
    if ccon > 0:
        trtConMat, trtConES = getDiffAggLoops(td, conLoops, cpu)
        sns.heatmap(trtConMat,
                    xticklabels=False,
                    yticklabels=False,
                    square=True,
                    ax=ax,
                    cmap=cmap,
                    linewidths=0.05,
                    linecolor="gray",
                    linestyle="--",
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={"shrink": 0.5})
        ax.set_title("%s specific loops\nES:%.3f" %
                     (ccon, np.mean(trtConES)),
                     fontsize=8)
    else:
        ax.set_title("No %s unique loops" % nb)

    ax = axs[1][0]
    if cover > 0:
        conOverMat, conOverES = getDiffAggLoops(cd, overlappedLoops, cpu)
        sns.heatmap(conOverMat,
                    xticklabels=False,
                    yticklabels=False,
                    square=True,
                    ax=ax,
                    cmap=cmap,
                    linewidths=0.05,
                    linecolor="gray",
                    linestyle="--",
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={"shrink": 0.5})
        ax.set_ylabel(nb, fontsize=10)
        ax.set_title("ES:%.3f" %np.mean(conOverES), fontsize=8)

    ax = axs[1][1]
    if ctrt > 0:
        conTrtMat, conTrtES = getDiffAggLoops(cd, trtLoops, cpu)
        sns.heatmap(conTrtMat,
                    xticklabels=False,
                    yticklabels=False,
                    square=True,
                    ax=ax,
                    cmap=cmap,
                    linewidths=0.05,
                    linecolor="gray",
                    linestyle="--",
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={"shrink": 0.5})
        ax.set_title("ES:%.3f" %np.mean(conTrtES),fontsize=8)

    ax = axs[1][2]
    if ccon > 0:
        conConMat, conConES = getDiffAggLoops(cd, conLoops, cpu)
        sns.heatmap(conConMat,
                    xticklabels=False,
                    yticklabels=False,
                    square=True,
                    ax=ax,
                    cmap=cmap,
                    linewidths=0.05,
                    linecolor="gray",
                    linestyle="--",
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={"shrink": 0.5})
        ax.set_title("ES:%.3f" % np.mean(conConES), fontsize=8)

    pylab.tight_layout()
    pylab.savefig(output + "_diffAggLoops.pdf")


def callDiffLoops(
        tl,
        cl,
        td,
        cd,
        output,
        cut=0,
        mcut=-1,
        cpu=1,
        pcut=1e-2,
        igp=False,
        noPCorr=False,
        fdrcut=0.05,
        juicebox=False,
        washU=False,
        customize=False,
        cacut=0.0,
        cmcut=0.0,
        vmin=None,
        vmax=None,
        cmap=None,
):
    """
    Call differentially enriched loops 
    @param tl: str, file of _loops.txt for treatment sample
    @param cl: str, file of _loops.txt for control sample
    @param td: str, directory generated by cLoops2 pre for treatment sample
    @param cd: str, directory generated by cLoops2 pre for control sample 
    @param output: str, prefix for output file 
    @param cut: int, distance cutoff for estimation of difference significance , >=cut
    @param mcut: int, distance cutoff for estimation of difference significance, <=mcut
    @param cpu: int, number of cpus used 
    @param pcut: float, p-value cutoffs after Bon correction
    @param igp: bool, whether to ignore p-value cutoff
    @param noPCorr: bool, whehter to perform Bon correction of p-values, default yes
    @param fdrcut: float, fdrcut for background to estimate Mcut and Acut
    @param customize: binary, if true, use user provided MA M cut and A cut
    @param cacut: float, if customize, used, A for MA plot
    @param cmcut: float, if customize, used, M for MA plot
    @param cmap: str, color map string option
    """
    #data name
    if td.endswith("/"):
        td = td[:-1]
    if cd.endswith("/"):
        cd = cd[:-1]
    tname = td.split("/")[-1]
    cname = cd.split("/")[-1]

    #read in loops
    tloops = parseTxt2Loops(tl)
    cloops = parseTxt2Loops(cl)

    #process meta information
    na = td.split("/")[-1]  #name of sample directory
    tmetaf = td + "/petMeta.json"
    tmeta = json.loads(open(tmetaf).read())
    nb = cd.split("/")[-1]
    cmetaf = cd + "/petMeta.json"
    cmeta = json.loads(open(cmetaf).read())

    #total PETs
    ta = tmeta["Unique PETs"]
    tb = cmeta["Unique PETs"]

    #chromosomes for testing
    keys = set(tmeta["data"]["cis"].keys()).intersection(
        set(cmeta["data"]["cis"].keys()))

    # step 1, merge the overlapped loops
    mloops = mergeLoops(tloops, na, cloops, nb)
    keys = list(keys.intersection(mloops.keys()))

    # step 2, quantify the loops in two conditions
    ds = Parallel(n_jobs=cpu, backend="multiprocessing")(delayed(quantDloops)(
        key,
        mloops[key],
        tmeta["data"]["cis"][key]["ixy"],
        cmeta["data"]["cis"][key]["ixy"],
        cut=cut,
        mcut=mcut,
    ) for key in keys)
    ts, cs = [], []
    dloops = {}
    for d in ds:
        if d is None:
            continue
        ts.extend(d[1])
        cs.extend(d[2])
        dloops[d[0]] = d[3]

    # step 3, estimate the fitting parameters, cutoffs based on MANorm
    sf, acut, mcut = getBgNorm(cs, ts, output, fdrcut=fdrcut)
    # check whether to use customized cutoffs
    if customize:
        acut = cacut 
        mcut = cmcut

    # step 4, estimate the difference significance
    ds = Parallel(n_jobs=cpu,
                  backend="multiprocessing")(delayed(estLoopDiffSig)(
                      key,
                      sf,
                      ta,
                      tb,
                      dloops[key],
                  ) for key in keys)
    dloops = []
    cs, ts, ts2 = [], [], []
    for d in ds:
        if d is None:
            continue
        dloops.extend(d[0])
        cs.extend(d[1])
        ts.extend(d[2])
        ts2.extend(d[3])

    #step 5, p-values Bonferroni correction and determine whether significant
    dloops = markDiffSig(dloops, acut, mcut, pcut=pcut,igp=igp,noPCorr=noPCorr)
    sigIndex = [i for i, loop in enumerate(dloops) if loop.significant > 0]

    # step 6, write the result
    dloops2txt(dloops, output + "_dloops.txt")

    # step 7, write the result as washU or juicebox
    tloops = [
        dloop for dloop in dloops
        if dloop.significant > 0 and dloop.scaled_fc > 0
    ]
    cloops = [
        dloop for dloop in dloops
        if dloop.significant > 0 and dloop.scaled_fc < 0
    ]
    dloops2txt( tloops, output + "_" + tname +"_specific_dloops.txt")
    dloops2txt( cloops, output + "_" + cname +"_specific_dloops.txt")
    comloops = [dloop for dloop in dloops if dloop.significant <1]
    if juicebox:
        dloops2juiceTxt(tloops, output + "_" + tname + "_loops_juicebox.txt")
        dloops2juiceTxt(cloops, output + "_" + cname + "_loops_juicebox.txt")
        dloops2juiceTxt(comloops, output + "_common_loops_juicebox.txt",significant=0)
    if washU:
        loops2washuTxt(tloops, output + "_" + tname + "_loops_legacyWashU.txt")
        loops2washuTxt(cloops, output + "_" + cname + "_loops_legacyWashU.txt")
        loops2washuTxt(comloops, output + "_common_loops_legacyWashU.txt",significant=0)
        dloops2NewWashuTxt(tloops,
                           output + "_" + tname + "_loops_newWashU.txt")
        dloops2NewWashuTxt(cloops,
                           output + "_" + cname + "_loops_newWashU.txt")
        dloops2NewWashuTxt(comloops, output + "_common_loops_newWashU.txt",significant=0)

    # step 8, show plot
    #ma plot
    plotDiffLoopsMA(sigIndex, cs, ts, ts2, tname, cname, mcut, acut, output)
    #volcano plot
    plotDiffLoopsVolcano(output + "_dloops.txt",
                         output,
                         tname,
                         cname,
                         fccut=mcut,
                         pcut=pcut)
    #plot aggregated differential loops
    plotDiffAggLoops(dloops, output, tl, cl, td, cd, cpu=cpu, norm=True,vmin=vmin,vmax=vmax,cmap=cmap)

