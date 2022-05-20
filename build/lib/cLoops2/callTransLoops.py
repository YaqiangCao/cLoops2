#!/usr/bin/env python3
#--coding:utf-8 --
"""
callTransLoops.py
"""

#sys
import os
import sys
import json
from glob import glob
from datetime import datetime
from collections import Counter

#3rd
import joblib
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn import linear_model
from scipy.stats import hypergeom, binom, poisson, combine_pvalues

#cLoops
from cLoops2.ds import Loop, XY
from cLoops2.io import parseIxy, ixy2pet, loops2juiceTxt, loops2washuTxt, ixy2pet, updateJson, loops2txt
from cLoops2.geo import checkLoopOverlap, combineLoops
from cLoops2.settings import *
from cLoops2.blockDBSCAN import blockDBSCAN as DBSCAN
from cLoops2.callCisLoops import getPerRegions, selSigLoops, estAnchorSig


def runTransDBSCANLoops(fixy, eps, minPts):
    """
    Run DBSCAN to detect interactions for one .ixy file.
    @param fixy: str, .ixy file name 
    @param eps: int, eps for DBSCAN
    @param minPts: int, minPts for DBSCAN
    """
    loops, loopReads = [], []
    key, mat = parseIxy(fixy, cut=0)
    #gave mat each PET a id
    mat2 = np.zeros((mat.shape[0],3))
    mat2[:,0] = np.arange(mat.shape[0])
    mat2[:,1] = mat[:,0]
    mat2[:,2] = mat[:,1]
    mat = mat2

    #data for interaction records, read for readId
    report = "%s \t Clustering %s and %s using eps %s, minPts %s\n" % (
        datetime.now(), key[0], key[1], eps, minPts)
    sys.stderr.write(report)
    db = DBSCAN(mat, eps, minPts)
    labels = pd.Series(db.labels)
    mat = pd.DataFrame(mat[:, 1:].astype("int"),
                       index=mat[:, 0],
                       columns=["X", "Y"])
    nlabels = set(labels.values)
    #collect clusters
    for label in nlabels:
        los = list(labels[labels == label].index)
        loopReads.extend(los)
        sub = mat.loc[los, :]
        if int(np.min(sub["X"])) == int(np.max(sub["X"])) or int(
                np.min(sub["Y"])) == int(np.max(sub["Y"])):
            continue
        #define loops
        loop = Loop()
        loop.rab = sub.shape[0]
        loop.chromX = key[0]
        loop.chromY = key[1]
        loop.x_start = int(np.min(sub["X"]))
        loop.x_end = int(np.max(sub["X"]))
        loop.x_center = (loop.x_start + loop.x_end) / 2
        loop.y_start = int(np.min(sub["Y"]))
        loop.y_end = int(np.max(sub["Y"]))
        loop.y_center = (loop.y_start + loop.y_end) / 2
        loop.distance = -1
        loop.cis = False
        loops.append(loop)
    report = "%s \t Clustering %s and %s finished. Estimated %s reads for %s candidate loops. \n" % (
        datetime.now(), key[0], key[1], len(loopReads), len(loops))
    sys.stderr.write(report)
    return "-".join(key), loops


#related
def parallelRunTransDBSCANLoops(meta, eps, minPts, cpu=1):
    """
    Paralle version of runCisDBSCANLoops
    @param meta: meta information parsed form petMeta.json
    @param eps: int, eps for DBSCAN
    @param minPts: int, minPts for DBSCAN
    """
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(runTransDBSCANLoops)(
        meta["data"]["trans"][key]["ixy"], eps, minPts)
                              for key in meta["data"]["trans"].keys())
    loops = {}
    for d in ds:
        if d is not None and len(d[1]) > 0:
            key, di = d[0], d[1]
            loops[key] = di
    return loops


def estLoopSig(key,
               loops,
               fixy,
               minPts=5,
               pseudo=1,
               peakPcut=1e-5,
               peakFccut=2,
               countDiffCut=10):
    """
    Estimate the loop statstical significance for one chromosomal.
    @param loops: list of Loop object
    @param fixy: cLoops2 pre generated .ixy file
    @param minPts: int, minPts
    """
    xy = ixy2pet(fixy, cut=0)
    N = xy.number
    print("%s \t Estimate significance for %s candidate interactions in %s." %
          (datetime.now(), len(loops), key))
    nloops = []
    for loop in tqdm(loops):
        ra, rb, rab = xy.queryLoop(loop.x_start, loop.x_end, loop.y_start,
                                   loop.y_end)
        ra, rb, rab = len(ra), len(rb), len(rab)
        if rab < minPts:
            continue
        loop.ra = ra
        loop.rb = rb
        loop.rab = rab
        #unbalanced anchor density, to avoid lines, unknow reason for lines, maybe stripes
        if ra / float(rb) > countDiffCut or rb / float(ra) > countDiffCut:
            continue
        if (loop.x_end -
                loop.x_start) / (loop.y_end - loop.y_start) > countDiffCut or (
                    loop.y_end - loop.y_start) / (loop.x_end -
                                                  loop.x_start) > countDiffCut:
            continue
        lowerra, lowerrb, lowerrab = xy.queryLoop(
            loop.x_start - (loop.x_end - loop.x_start), loop.x_start,
            loop.y_start - (loop.y_end - loop.y_start), loop.y_start)  #p2ll
        loop.P2LL = float(rab) / max(len(lowerrab), pseudo)
        px, esx = estAnchorSig(xy, loop.x_start, loop.x_end)
        py, esy = estAnchorSig(xy, loop.y_start, loop.y_end)
        loop.x_peak_poisson_p_value = px
        loop.x_peak_es = esx
        loop.y_peak_poisson_p_value = py
        loop.y_peak_es = esy
        #hypergeometric p-value
        hyp = max([1e-300, hypergeom.sf(rab - 1.0, N, ra, rb)])
        #start caculate the permutated background
        rabs, nbps = [], []
        nas, nbs = getPerRegions(loop, xy)
        for na in nas:
            nac = float(len(na))
            for nb in nbs:
                nbc = float(len(nb))
                nrab = float(len(na.intersection(nb)))
                #collect the value for poisson test
                rabs.append(nrab)
                #collect the possibility for following binomial test
                if nac > 0 and nbc > 0:
                    den = nrab / (nac * nbc)
                    nbps.append(den)
                else:
                    nbps.append(0.0)
        rabs, nbps = np.array(rabs), np.array(nbps)
        mrabs = float(np.mean(rabs))
        mbps = np.mean(nbps)
        #local fdr
        fdr = len(rabs[rabs > rab]) / float(len(rabs))
        #enrichment score
        es = rab / max(mrabs, pseudo)
        #simple possion test
        pop = max([1e-300, poisson.sf(rab - 1.0, mrabs)])
        #simple binomial test
        nbp = max([
            1e-300, binom.sf(rab - 1.0, ra * rb, mbps)
        ])  #the p-value is quit similar to that of cLoops 1 binomial test
        loop.FDR = fdr
        loop.ES = es
        loop.density = float(
            loop.rab) / (loop.x_end - loop.x_start + loop.y_end -
                         loop.y_start) / N * 10.0**9
        loop.hypergeometric_p_value = hyp
        loop.poisson_p_value = pop
        loop.binomial_p_value = nbp
        nloops.append(loop)
        #print(ra,rb,rab,mrabs,es,fdr,hyp,pop,nbp,n,nbp2)
    return key, nloops


def markSigLoops(key, loops):
    """
    Mark the significance of different loops.
    """
    sig = lambda x: True if x.binomial_p_value <= 1e-10 and x.FDR <= 0.05 and loop.ES >= 2 else False
    for loop in loops:
        if sig(loop):
            loop.significant = 1
        else:
            loop.significant = 0
    return key, loops


def callTransLoops(
        predir,
        fout,
        logger,
        eps=[2000, 5000],
        minPts=[5, 10],
        cpu=1,
        filter=False,
        washU=False,
        juicebox=False,
):
    """
    Call inter-chromosomal loops parallel.
    @param metaf: str, petMeta.json file for calling peaks
    @param eps: list
    @param minPts: list
    """

    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    ## step 1 find the candidate loops by running multiple times of clustering
    loops = {}  #candidate loops
    for ep in eps:
        for minPt in minPts:
            loops_2 = parallelRunTransDBSCANLoops(meta, ep, minPt, cpu=cpu)
            loops = combineLoops(loops, loops_2)
    ## step 2 determine the statstical significance of candidate loops
    logger.info("Estimating loop statstical significance.")
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(estLoopSig)(
        key,
        loops[key],
        meta["data"]["trans"][key]["ixy"],
        minPts=max(minPts),
    ) for key in loops.keys())
    nds = {}
    for d in ds:
        nds[d[0]] = d[1]

    #mark the significant loops
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(markSigLoops)(key, nds[key])
                              for key in nds.keys())
    nds = {}
    for d in ds:
        nds[d[0]] = d[1]

    ## step 4 for the overlapped loops, output the most significant one
    logger.info("Selecting the most significant loops of overlapped ones. ")
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(selSigLoops)(key, nds[key])
                              for key in nds.keys())
    nds = {}
    for d in ds:
        nds[d[0]] = d[1]
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(selSigLoops)(key, nds[key])
                              for key in nds.keys())
    nds = {}
    for d in ds:
        nds[d[0]] = d[1]
    loops = []
    for d in ds:
        loops.extend(d[1])
    ## step 5 output
    logger.info("Output %s loops to %s_loops.txt" % (len(loops), fout))
    loops2txt(loops, fout + "_trans_loops.txt")
    if washU:
        loops2washuTxt(loops, fout + "_trans_loops_washU.txt")
    if juicebox:
        loops2juiceTxt(loops, fout + "_trans_loops_juicebox.txt")
