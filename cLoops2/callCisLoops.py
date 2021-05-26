#!/usr/bin/env python3
#--coding:utf-8 --
"""
callCisLoops.py
2019-09-10: basically finished.
2019-09-26: due to low speed of XY object, now change the object way to function way. Do not kown why is so slow. Still same slow, change back to object method. The slow is due to blockDBSCAN called too many and too broad loops.
2019-09-29: updated binomial test method, remove fitting process and cut option from estLoopSig
2020-01-20: fine tune some details
2020-01-23: fine tune. 1) for trac-looping like data, binomial < 1e-1 is enough; 2) max_cut can speed up a lot for significance test and only few loops will lost. 3) blockDBSCAN is less sensitive to minPts. 4) filter loops with estimated distance cutoff, can be run within that loop. Out side with -max_cut may remove too many
2020-02-06: fine tune functions, getLoopNearbyPETs added.
2020-02-09: fine tune the anchor peak estimation, to loose mode
2020-02-12: change the nearby permutated PETs to median from mean, could be much sensitive to the result. Additionlly, will be less sensitive to eps, seems much better.
2020-02-13: both for TrAC-looping and HiC, blockDBSCAN is much snesitive and faster than cDBSCAN, so maybe therefore no more test 
2020-02-14: change HiC p2llcut to 1, much sensitive. 
2020-02-15: P2LL quite useless in cLoops2, no more trying. Finally binomial p value can control how many significant loops for HiC
2020-03-04: replace print to logger
2020-03-09: update density with library size, not the PETs number in that chromosome, more stable, not affect by estimated cutoffs
2020-03-11: remove the pseudo for estimate p-values of loops, for Trac-looping, it could at least >= 6 if pseudo =1 for poisson p < 1e-6, make the setting of minPts meanless
2020-11-22: using cDBSCAN2 for Hi-C data
2020-11-25: observed from Hi-C data, for overlapped loops, higher enrichment score,better
2021-03-23: change HIC P2LLcut to 1 and binomial p-value cut to 1e-3 as using cDBSCAN2; previouse cutoffs for HIC P2LLcut >=2 binomial p<=1e-5
2021-05-20: try to speed up permutation background query speed; tested with K562 Hi-TrAC chr21, 5 fold speed up.
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
from scipy.stats import hypergeom, binom, poisson

#cLoops
from cLoops2.settings import *
from cLoops2.ds import Loop, XY
from cLoops2.est import estIntraCut
from cLoops2.plot import plotIntraCut
#from cLoops2.blockDBSCAN import blockDBSCAN as DBSCAN
from cLoops2.geo import checkLoopOverlap, combineLoops
from cLoops2.io import parseIxy, ixy2pet, loops2juiceTxt, loops2washuTxt, updateJson, loops2txt, loops2ucscTxt,loops2NewWashuTxt

#gloabl settings 
logger = None
DBSCAN = None

def runCisDBSCANLoops(fixy, eps, minPts, cut=0,mcut=-1):
    """
    Run DBSCAN to detect interactions for one .ixy file.
    @param fixy: str, .ixy file name 
    @param eps: int, eps for DBSCAN
    @param minPts: int, minPts for DBSCAN
    """
    loops, loopReads, peakReads, distalDistance, closeDistance = [],[], [],[], []
    key, mat = parseIxy(fixy, cut=cut,mcut=mcut)
    mat2 = np.zeros((mat.shape[0],3))
    mat2[:,0] = range(mat.shape[0])
    mat2[:,1] = mat[:,0]
    mat2[:,2] = mat[:,1]
    mat = mat2
    mat = mat.astype("int")
    if key[0] != key[1]:
        return None
    if cut > 0:
        d = mat[:, 2] - mat[:, 1]
        p = np.where(d >= cut)[0]
        mat = mat[p, :]
        closeDistance.extend(list(d[d < cut]))
    if len(mat) == 0:
        report = "No PETs found in %s, maybe due to cut > %" % (fixy, cut)
        #print(report)
        logger.info(report)
        return None  #no data to call loops
    #data for interaction records, read for readId
    #report = "%s \t Clustering %s and %s using eps %s, minPts %s,pre-set distance cutoff > %s\n" % ( datetime.now(), key[0], key[1], eps, minPts, cut)
    #sys.stderr.write(report)
    report = "Clustering %s and %s using eps %s, minPts %s,pre-set distance cutoff > %s" % (key[0], key[1], eps, minPts, cut)
    logger.info(report)
    db = DBSCAN(mat, eps, minPts)
    labels = pd.Series(db.labels)
    mat = pd.DataFrame(mat[:, 1:].astype("int"),
                       index=mat[:, 0],
                       columns=["X", "Y"])
    nlabels = set(labels.values)
    #collect clusters
    for label in nlabels:
        los = list(labels[labels == label].index)
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
        loop.cis = True
        loop.distance = abs(loop.y_center - loop.x_center)
        if loop.x_end < loop.y_start:  #true candidate loops
            loops.append(loop)
            loopReads.extend(los)
        else:  #maybe peaks
            peakReads.extend(los)
    report = "Clustering %s and %s finished. Estimated %s self-ligation reads and %s inter-ligation reads, %s candidate loops." % (key[0], key[1], len(peakReads), len(loopReads), len(loops))
    logger.info(report)
    if len(loopReads) > 0:
        distalDistance = list(mat.loc[loopReads, "Y"] -
                              mat.loc[loopReads, "X"])
    if len(peakReads) > 0:
        closeDistance.extend(
            list(mat.loc[peakReads, "Y"] - mat.loc[peakReads, "X"]))
    return "-".join(key), loops, distalDistance, closeDistance


def parallelRunCisDBSCANLoops(meta, eps, minPts, cut=0,mcut=-1,cpu=1):
    """
    Paralle version of runCisDBSCANLoops
    @param meta: meta information parsed form petMeta.json
    @param eps: int, eps for DBSCAN
    @param minPts: int, minPts for DBSCAN
    """
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(runCisDBSCANLoops)(
        meta["data"]["cis"][key]["ixy"], eps, minPts, cut=cut,mcut=mcut)
                              for key in meta["data"]["cis"].keys())
    loops, dis, dss = {}, [], []
    for d in ds:
        if d is not None and len(d[1]) > 0:
            key, di, ddis, ddss = d[0], d[1], d[2], d[3]
            loops[key] = di
            dis.extend(ddis)
            dss.extend(ddss)
    return loops, dis, dss


def filterLoopsByDis(loops, cut):
    """
    Filter candidate loops by distance cutoffs
    """
    for key in loops:
        nr = []
        for loop in loops[key]:
            if loop.distance > cut:
                nr.append(loop)
        loops[key] = nr
    return loops


def getPerRegions(loop, xy, win=5):
    """
    Get the nearby regions for interacting two locus, win as how many nearby, 6 is enough for interacting more than 100 regions to estimate FDR and others. The mean distance of all the permutated regions is the same to that between iva and ivb.
    @param loop: cLoops2:ds:Loop 
    @param xy: cLoops2:ds:XY
    """
    ca = loop.x_center
    cb = loop.y_center
    sa = (loop.x_end - loop.x_start) / 2
    sb = (loop.y_end - loop.y_start) / 2
    nas, nbs = [], []
    step = (sa + sb) / 2
    #the permutated region all PET ids
    start = min([ ca-win*step-sa, cb-win*step-sb ])
    end = max([ca+win*step+sa,cb+win*step+sb])
    ps = list(xy.queryPeak( start, end))
    nmat = xy.mat[ps,]
    nxy = XY(nmat[:,0],nmat[:,1])
    # the PET id in the permutated regions
    for i in range(0 - win, win + 1):
        if i == 0:
            continue
        niva = [max([0, ca + i * step - sa]), max([0, ca + i * step + sa])]
        nivb = [max([0, cb + i * step - sb]), max([0, cb + i * step + sb])]
        nas.append(nxy.queryPeak(niva[0], niva[1]))
        nbs.append(nxy.queryPeak(nivb[0], nivb[1]))
    return nas, nbs


def getLoopNearbyPETs(loop, xy, win=5):
    """
    Get the target loop nearby PETs
    """
    nas, nbs = getPerRegions(loop, xy, win=win)
    rabs, nbps = [], []
    for na in nas:
        nac = float(len(na))
        for nb in nbs:
            nbc = float(len(nb))
            nrab = float(len(na.intersection(nb)))
            #collect the value for poisson test and binomial test
            if nrab > 0:
                rabs.append(nrab)
                den = nrab / (nac * nbc)
                nbps.append(den)
            #skip zeros will be more strong test, for Trac-looping will be no significant loops
            #need to test for Hi-C if remove 0 will be much better
            else:
                rabs.append(0)
                nbps.append(0.0)
    return np.array(rabs), np.array(nbps)


def estAnchorSig(xy, left, right, ext=5):
    """
    Estimate the anchor significance as peak, using the similar idea of callPeaks.
    """
    rpb = float(xy.number) / (np.max(xy.ys) - np.min(xy.xs))
    m = (left + right) / 2
    length = right - left
    #if using queryPeakBoth, no significant anchors, do not try to use queryPeakBoth again
    count = len(xy.queryPeak(left, right))
    cs = []
    #local test
    start = max(0, m - ext * length)
    end = m + ext * length
    r = (len(xy.queryPeak(start, end)) -
         count) / ext / 2  #mean value of nearby
    cs.append(r)
    #global test
    cs.extend([rpb * length, 1])  #1 is used as pesudo count
    c = float(max(cs))
    es = count / c
    p = max([1e-300, poisson.sf(count - 1.0, c)])
    return p, es


def estLoopSig(
        key,
        loops,
        fixy,
        tot,
        minPts=5,
        pseudo=1,
        cut=0,
        mcut=-1,
        peakPcut=1e-5,
        win=5,
        countDiffCut=20,
        hic=False,
):
    """
    Estimate the loop statstical significance for one chromosomal.
    @param loops: list of Loop object
    @param fixy: cLoops2 pre generated .ixy file
    @param hic: bool, if True, will skip anchor examazaiton and carry P2LL examazation
    """
    if hic:
        p2llcut = 1
    else:
        p2llcut = 1 
    xy = ixy2pet(fixy, cut=cut,mcut=mcut)
    N = xy.number
    logger.info( "Estimate significance for %s candidate interactions in %s with %s PETs distance > =%s and <=%s,requiring minPts >=%s." % (len(loops), key, N, cut, mcut,minPts))
    nloops = []
    for loop in tqdm(loops):
        #filtering unbalanced anchor size
        if (loop.x_end -
                loop.x_start) / (loop.y_end - loop.y_start) > countDiffCut or (
                    loop.y_end - loop.y_start) / (loop.x_end -
                                                  loop.x_start) > countDiffCut:
            continue
        ra, rb, rab = xy.queryLoop(loop.x_start, loop.x_end, loop.y_start,
                                   loop.y_end)
        ra, rb, rab = len(ra), len(rb), len(rab)
        if rab < minPts:
            continue
        #unbalanced anchor density, to avoid lines, unknow reason for lines, maybe stripes
        if ra / float(rb) > countDiffCut or rb / float(ra) > countDiffCut:
            continue
        loop.ra = ra
        loop.rb = rb
        loop.rab = rab
        #P2LL
        lowerra, lowerrb, lowerrab = xy.queryLoop(
            loop.x_end,  (loop.x_end - loop.x_start) + loop.x_end,
            loop.y_start - (loop.y_end - loop.y_start),
            loop.y_start)  #p2ll, seems useless
        loop.P2LL = float(rab) / max(len(lowerrab), pseudo)
        if hic and loop.P2LL < p2llcut:
            continue
        #hypergeometric p-value, if the simple hy test can not pass, no need for furthur test
        hyp = max([1e-300, hypergeom.sf(rab - 1.0, N, ra, rb)])
        if hyp > 1e-2:
            continue
        #start caculate the permutated background
        rabs, nbps = getLoopNearbyPETs(loop, xy, win)
        mrabs = float(np.median(rabs))
        mbps = np.median(nbps)
        #local fdr
        if len(rabs) > 0:
            fdr = len(rabs[rabs > rab]) / float(len(rabs))
        else:
            fdr = 0.0
        if mrabs >= rab or fdr >= 0.1:  #hope to speed up
            continue
        #enrichment score
        es = rab / max(mrabs,pseudo) #pseudo only used to avoid inf
        if es < 1: #hope to speed up 
            continue
        #simple possion test
        pop = max([1e-300, poisson.sf(rab - 1.0, mrabs)])
        #simple binomial test
        nbp = max([
            1e-300, binom.sf(rab - 1.0, ra * rb - rab, mbps)
        ])  #the p-value is quit similar to that of cLoops 1 binomial test
        #nbp = max([1e-300, binom.sf(rab - 1.0, N - rab, mbps * ra * rb / N)])  #cLoops 1 binomial test
        loop.FDR = fdr
        loop.ES = es
        loop.density = float(
            loop.rab) / (loop.x_end - loop.x_start + loop.y_end -
                         loop.y_start) / tot * 10.0**9
        loop.hypergeometric_p_value = hyp
        loop.poisson_p_value = pop
        loop.binomial_p_value = nbp

        #make sure the anchor are significant
        px, esx = estAnchorSig(xy, loop.x_start, loop.x_end)
        py, esy = estAnchorSig(xy, loop.y_start, loop.y_end)
        if hic == False and not (px < peakPcut and py < peakPcut):
            continue
        loop.x_peak_poisson_p_value = px
        loop.x_peak_es = esx
        loop.y_peak_poisson_p_value = py
        loop.y_peak_es = esy
        nloops.append(loop)
    return key, nloops


def markSigLoops(key, loops, hic=False):
    """
    Mark the significance of different loops.
    """
    sig = lambda x: True if x.FDR <= 0.05 and x.ES >= 2 and x.hypergeometric_p_value <= 1e-5 and x.poisson_p_value <= 1e-5 else False
    for loop in loops:
        if sig(loop):
            if hic:
                if loop.binomial_p_value < 1e-3:
                    loop.significant = 1
                else:
                    loop.significant = 0
            else:
                if loop.binomial_p_value < 1e-1:
                    loop.significant = 1
                else:
                    loop.significant = 0
        else:
            loop.significant = 0
    return key, loops


def selSigLoops(key, loops):
    """
    Remove overlapped called loops, keep the more significant one for multiple eps result. 
    """
    #only consider the significant loops to reduce search space
    loops = [loop for loop in loops if loop.significant > 0]
    #find the overlapped loops
    nloops = []
    skips = set()
    for i in range(len(loops)):
        if i in skips:
            continue
        n = [loops[i]]
        for j in range(i + 1, len(loops)):
            for p in n:
                if checkLoopOverlap(p, loops[j]):
                    n.append(loops[j])
                    skips.add(j)
                    break
        nloops.append(n)
    #get the most significant loops of the overlapped ones according to enrichment score
    nnloops = []
    for n in nloops:
        if len(n) == 1:
            nnloops.append(n[0])
        else:
            for i in range(len(n) - 1):
                for j in range(i + 1, len(n)):
                    #if n[i].binomial_p_value > n[j].binomial_p_value:
                    if n[i].ES < n[j].ES: #these options actually does not matter a lot, observed from Hi-C
                    #if n[i].density < n[j].density:
                        n[i], n[j] = n[j], n[i]
            nnloops.append(n[0])
    #search again, in case of lost
    for loopa in loops:
        flag = 0
        for i, loopb in enumerate(nnloops):
            if checkLoopOverlap(loopa, loopb):
                flag = 1
                break
        if flag == 0:
            nnloops.append(loopa)
    return key, nnloops


def getAllAnchors(loops, margin=1):
    """
    Get the genomic set of all anchors.
    """
    cov = set()
    for loop in loops:
        cov.update(range(loop.x_start, loop.x_end + 1))
        cov.update(range(loop.y_start, loop.y_end + 1))
    cov = list(cov)
    cov.sort()
    anchors = []
    i = 0
    while i < len(cov) - 1:
        for j in range(i + 1, len(cov)):
            if cov[j] - cov[j - 1] > margin:
                break
            else:
                continue
        start = cov[i]
        end = cov[j - 1]
        anchors.append([start, end])
        i = j  #update search start
    return anchors


def filterPETs(key, predir, fixy, loops, margin=1):
    """
    Filter PETs, only keep those located at loop anchor regions.
    """
    #print("%s\t Filtering PETs of %s with %s loops." % (datetime.now(), key, len(loops)))
    logger.info("Filtering PETs of %s with %s loops." % (key, len(loops)))
    anchors = getAllAnchors(loops, margin=margin)
    key2, mat = parseIxy(fixy)
    xy = XY(mat[:,0],mat[:,1]) 
    rs = set()
    for iv in anchors:
        r = xy.queryPeak(iv[0], iv[1])
        rs.update(r)
    rs = list(rs)
    if len(rs) == 0:
        return
    mat = mat[rs, ]
    foixy = predir + "/" + "-".join(key2) + ".ixy"
    joblib.dump(mat, foixy)


def callCisLoops(
        predir,
        fout,
        log,
        eps=[2000, 5000],
        minPts=[5, 10],
        cpu=1,
        cut=0,
        mcut=-1,
        plot=False,
        max_cut=False,
        hic=False,
        filter=False,
        ucsc=False,
        juicebox=False,
        washU=False,
        emPair=False,
):
    """
    Call intra-chromosomal loops parallel.
    @param metaf: str, petMeta.json file for calling peaks
    @param eps: list
    @param minPts: list
    @param empair: bool, if true, pair run eps and minPts, 
    """
    global logger
    logger = log
    global DBSCAN 
    if hic:
        from cLoops2.cDBSCAN2 import cDBSCAN as DBSCAN
        logger.info("-hic option selected, cDBSCAN2 is used instead of blockDBSCAN.")
    else:
        from cLoops2.blockDBSCAN import blockDBSCAN as DBSCAN
    if emPair and len(eps) != len(minPts):
        logger.info("-emPair option selected, number of eps not equal to minPts, return.")
        return 
    ##step 0 prepare data and check directories
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    tot = meta["Unique PETs"]
    if filter:
        logger.info(
            "-filter option chosed, will filter raw PETs based on called loops, for any PET that any end overlaps loop anchors will be kept. "
        )
        fdir = fout + "_filtered"
        if not os.path.exists(fdir):
            os.mkdir(fdir)
        elif len(os.listdir(fdir)) > 0:
            r = "working directory %s exists and not empty." % fdir
            logger.error(r)
            return

    ## step 1 find the candidate loops by running multiple times of clustering
    loops = {}  #candidate loops
    #distance of classified inter-ligation PETs, self-ligaiton PETs.
    dis, dss = [], []
    cuts = [
        cut,
    ]
    if emPair:
        for ep,minPt in zip(eps,minPts):
            loops_2, dis_2, dss_2 = parallelRunCisDBSCANLoops(
                    meta,
                    ep,
                    minPt,
                    cut=cut,
                    mcut=mcut,
                    cpu=cpu,
                )
            if len(dis_2) == 0:
                logger.error(
                    "ERROR: no inter-ligation PETs detected for eps %s minPts %s,can't model the distance cutoff,continue anyway"
                    % (ep, minPt))
                continue
            if not (len(dis_2) == 0 or len(dss_2) == 0):
                cut_2 = estIntraCut(np.array(dis_2), np.array(dss_2))
                if plot:
                    plotIntraCut(dis_2,
                                 dss_2,
                                 cut_2,
                                 prefix=fout + "_eps%s_minPts%s_disCutoff" %
                                 (ep, minPt))
                logger.info(
                    "Estimated inter-ligation and self-ligation distance cutoff > %s for eps=%s,minPts=%s"
                    % (cut_2, ep, minPt))
            if len(dss_2) == 0:
                logger.info(
                    "No self-ligation PETs found, using cutoff > %s for eps=%s,minPts=%s"
                    % (cut, ep, minPt))
                cut_2 = cut
            loops_2 = filterLoopsByDis(loops_2, cut_2)
            loops = combineLoops(loops, loops_2)
            cuts.append(cut)
            cut = cut
    else:
        for ep in eps:
            for minPt in minPts:
                loops_2, dis_2, dss_2 = parallelRunCisDBSCANLoops(
                    meta,
                    ep,
                    minPt,
                    cut=cut,
                    mcut=mcut,
                    cpu=cpu,
                )
                if len(dis_2) == 0:
                    logger.error(
                        "ERROR: no inter-ligation PETs detected for eps %s minPts %s,can't model the distance cutoff,continue anyway"
                        % (ep, minPt))
                    continue
                if not (len(dis_2) == 0 or len(dss_2) == 0):
                    cut_2 = estIntraCut(np.array(dis_2), np.array(dss_2))
                    if plot:
                        plotIntraCut(dis_2,
                                     dss_2,
                                     cut_2,
                                     prefix=fout + "_eps%s_minPts%s_disCutoff" %
                                     (ep, minPt))
                    logger.info(
                        "Estimated inter-ligation and self-ligation distance cutoff > %s for eps=%s,minPts=%s"
                        % (cut_2, ep, minPt))
                if len(dss_2) == 0:
                    logger.info(
                        "No self-ligation PETs found, using cutoff > %s for eps=%s,minPts=%s"
                        % (cut, ep, minPt))
                    cut_2 = cut
                loops_2 = filterLoopsByDis(loops_2, cut_2)
                loops = combineLoops(loops, loops_2)
                cuts.append(cut_2)
                cut = cut_2

        #distance cutoff for estimation of loop significance
    #cuts = [c for c in cuts if c > 0]
    ncuts = [c for c in cuts if c > cuts[0]]
    ncuts.append( cuts[0] )
    cuts = ncuts
    if max_cut:
        cut = np.max(cuts)
    else:
        cut = np.min(cuts)

    ## step 2 determine the statstical significance of candidate loops
    logger.info("Estimating loop statstical significance.")
    if emPair:
        mm = min(minPts)
    else:
        mm = max(minPts)
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(
        delayed(estLoopSig)(
            key,
            loops[key],
            meta["data"]["cis"][key]["ixy"],
            tot,
            #minPts=max(minPts),
            minPts=mm,
            #cut= 0,  #if using estimated cut, will generate just a little few loops than cut=0, but will increase a lot speed
            cut=cut,
            mcut=mcut,
            hic=hic) for key in loops.keys())
    nds = {}
    for d in ds:
        nds[d[0]] = d[1]

    #mark the significant loops
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(markSigLoops)(key, nds[key], hic=hic)
                              for key in nds.keys())
    nds = {}
    for d in ds:
        nds[d[0]] = d[1]

    ## step 3 for the overlapped loops, output the most significant one
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

    ## step 4 output
    logger.info("Output %s loops to %s_loops.txt" % (len(loops), fout))
    loops2txt(loops, fout + "_loops.txt")
    if ucsc:
        loops2ucscTxt(loops, fout + "_loops_ucsc.interact")
    if juicebox:
        loops2juiceTxt(loops, fout + "_loops_juicebox.txt")
    if washU:
        loops2washuTxt(loops, fout + "_loops_legacyWashU.txt")
        loops2NewWashuTxt(loops, fout + "_loops_newWashU.txt")

    ## step 5 filtering PETs according to called loops
    if filter:
        Parallel(n_jobs=cpu,backend="multiprocessing")(
            delayed(filterPETs)(key,
                                fdir,
                                meta["data"]["cis"][key]["ixy"],
                                nds[key],
                                margin=max(eps)) for key in nds.keys())
        ixyfs = glob(fdir + "/*.ixy")
        tot = 0
        for f in ixyfs:
            key, mat = parseIxy(f)
            tot += mat.shape[0]
        nmetaf = fdir + "/petMeta.json"
        with open(nmetaf, "w") as fo:
            json.dump({"Unique PETs": tot}, fo)
        updateJson(ixyfs, nmetaf)
