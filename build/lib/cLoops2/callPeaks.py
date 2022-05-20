#!/usr/bin/env python3
#--coding:utf-8 --
"""
callPeaks.py
2019-08-27: updated as select the most significant peaks for overlapped ones; also stich together close peaks.
2019-09-10: basically finished.
2020-01-20: fine tune, also change cDBSCAN to blockDBSCAN
2020-01-25: fine tune enrichment score, using real control RPKM, only output significant peaks 
2020-01-25: sequencing depth ratio as scaling factor, enrichment compared to nearby cutoff added
2020-02-12: fine tune nearby regions, change from mean to median, indeed more sensitive, maybe min(mean,median) will be better to handle nearby there is a big peak. Fix the sensitive problem, by requring the queryPeakBoth in each nearby same size windows.
2020-02-13: the algorithm nearly perfect
2020-03-04: change print to logger.info for well organized log info, with setting logger as global var
2020-03-11: update fine tune for control local background
2020-03-20: change the linear fitting using nearby reads,if using the condidate peak region, then the assumption is all the reads in peaks are potential noise.
2020-04-01: fine tune compared to bg, after compared to bg, then do select signficiant peaks for overlapped.
2020-04-06: change the compare to nearby regions from meian to mean, more conserved.
2020-04-07: p-values correction to q-values added; harmonic mean for p_values added:https://www.pnas.org/content/116/4/1195
2020-04-08: adjusted p-values as final peaks selecting metric
2020-04-23: sensitive mode added. 
2020-07-29: -cut and -mcut integrated. -split added, only using single-end as bed for calling peaks, sometimes usefule for Trac-looping or HiChIP
2020-12-17: refined re-search of significant peaks.
2021-04-01: adding summit
"""

#sys
import json
import sys
from copy import deepcopy
from datetime import datetime

#3rd
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import poisson
from joblib import Parallel, delayed

#cLoops
from cLoops2.est import estSf
from cLoops2.ds import Peak, XY
from cLoops2.geo import stichPeaks, checkPeakOverlap
from cLoops2.io import parseIxy, ixy2pet, peaks2txt, peaks2bed
from cLoops2.blockDBSCAN import blockDBSCAN as DBSCAN
from cLoops2.cmat import get1DSig

#gloabl settings
logger = None

def getSplitMat(mat,splitExt):
    """
    Get the paired-end splited single end reads.
    """
    nmat = []
    for t in mat:
        nmat.append( [t[0]-splitExt,t[0]+splitExt] )
        nmat.append( [t[1]-splitExt,t[1]+splitExt] )
    mat = np.array(nmat)
    return mat
 

def runDBSCANPeaks(fixy, eps, minPts, cut=0,mcut=-1,split=False,splitExt=50):
    """
    Run DBSCAN to detect interactions for one .ixy file.
    @param f: str, .ixy file name 
    @param eps: int, eps for DBSCAN
    @param minPts: int, minPts for DBSCAN
    @return key: tuple,(chrA, chrB)
    @return peaks: list, peaks object
    """
    peaks, readS = [], 0
    #limit the PETs distace < 1000 for calling peaks from interaction data
    key, mat = parseIxy(fixy, cut=cut,mcut=mcut)
    if key[0] != key[1]:
        return None
    if len(mat) == 0:
        report = "No PETs found in %s." % (fixy)
        logger.info(report)
        return None  #no data to call peaks
    #if split, get splited data
    if split:
        mat = getSplitMat( mat, splitExt)
    #gave mat each PET a id
    mat2 = np.zeros((mat.shape[0], 3))
    mat2[:, 0] = range(mat.shape[0])
    mat2[:, 1] = mat[:, 0]
    mat2[:, 2] = mat[:, 1]
    mat = mat2
    mat = mat.astype("int")
    #data for interaction records, read for readId
    #report = "%s \t Clustering %s to find candidate peaks using eps as %s, minPts as %s.\n" % ( datetime.now(), key[0], eps, minPts)
    #sys.stderr.write(report)
    report = "Clustering %s to find candidate peaks using eps as %s, minPts as %s." % (
        key[0], eps, minPts)
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
        x_start = int(np.min(sub["X"]))
        x_end = int(np.max(sub["X"]))
        if x_start == x_end:
            x_end += 1
        y_start = int(np.min(sub["Y"]))
        y_end = int(np.max(sub["Y"]))
        if y_start == y_end:
            y_end += 1
        #if x_end >= y_start and len(los) >= minPts:  #overlaped anchors
        if x_end + eps >= y_start and len(los) >= minPts:
            #overlaped anchors are peaks
            #define peaks
            peak = Peak()
            peak.chrom = key[0]
            peak.start = x_start
            peak.end = y_end
            peak.length = peak.end - peak.start + 1
            peaks.append(peak)
            readS += len(los)
    #report = "%s \t Clustering %s finished. Estimated %s PETs in candidate %s peaks. \n" % ( datetime.now(), key[0], readS, len(peaks))
    #sys.stderr.write(report)
    report = "Clustering %s finished. Estimated %s PETs in candidate %s peaks." % (
        key[0], readS, len(peaks))
    logger.info(report)
    report = "Stiching %s candidate close %s peaks." % (key[0], len(peaks))
    logger.info(report)
    #peaks = stichPeaks(peaks, margin=eps / 2)
    peaks = stichPeaks(peaks, margin=eps)
    logger.info("Merged %s %s peaks" % (key[0], len(peaks)))
    return peaks


def removeSamePeaks(peaks):
    """
    Remove same peaks for different combined  eps and minPts clutering.
    """
    keys = set()
    uniquePeaks = []
    for peak in peaks:
        key = (peak.chrom, peak.start, peak.end)
        if key not in keys:
            keys.add(key)
            uniquePeaks.append(peak)
        else:
            continue
    return uniquePeaks


def sortPeaksByLength(peaks):
    """
    Sort peaks from longer to shoter.
    """
    #sort peaks according to size, from small to big
    peaks = sorted( peaks, key=lambda peak: peak.length )
    return peaks


def getPeakNearbyPETs(start, length, xy, win=5):
    """
    Get peak nearby background PETs numbers.
    """
    #end = peak.end + 1 + ext * peak.length
    cs = []
    for i in range(0 - win, win + 1):
        #if include itself, will be more conserved, but not good for clean data such as Cut&Tag
        if i == 0:
            continue
        s = max(0, start - 1 + i * length)
        e = max(0, start - 1 + (i + 1) * length)
        cs.append(len(xy.queryPeakBoth(s, e)))
    cs = np.array(cs)
    #cs = cs[cs > 0]  #with remove 0, only affect a little
    return cs


def findSigPeaks(
        fixy,
        totalPets,
        eps=[100, 200],
        minPts=[3, 5],
        fixybg=None,
        totalControlPets=None,
        pcut=1e-2,
        exts=[5, 10],
        lencut=200,
        escut=2.0,
        ESvsControlCut=2.0,
        pseudo=1.0,
        sen=False,
        cut=0,
        mcut=-1,
        split=False,
        splitExt=50,
):
    """
    Determine the significance of peaks based on poisson test for one chromosome.
    @param fixy: str, .ixy file
    @param totalPts: int, total PETs number
    @param eps: list of int, eps for DBSCAN
    @param minPts: list of int, minPts for DBSCAN
    @param fixybg: str, .ixy file of the control
    @param pcut: float, p-value cutoffs
    @param exts: list, extesion fold from the center of peaks to get significant test p-value, do not change if from outside, only used inside to debug.
    @param lencut: int, length cutoff for peaks.
    @pram escut: float, enrichment score cutoffs for peaks.
    @param ESvsControlCut: float, cutoff for enrichment score over control region
    @param pseudo: int, pseudo count to model background noise 
    @param sen: bool, whether use snesitive model to call peaks
    @param cut: int, distance >=cut PETs will be used
    @param mcut: int, distance <=mcut PETs will be used
    @param split: bool, whether to split paired-end tags as single-end reads
    @param splitExt: int, extension of the splited single-end read
    """
    peaks = []
    for ep in eps:
        for minPt in minPts:
            speaks = runDBSCANPeaks(fixy, ep, minPt,cut=cut,mcut=mcut,split=split,splitExt=splitExt)
            if speaks is None:
                continue
            peaks.extend(speaks)
    if len(peaks) == 0:
        return None
    #remove identical peaks
    peaks = removeSamePeaks(peaks)

    key, mat = parseIxy(fixy, cut=cut, mcut=mcut)
    if split:
        mat = getSplitMat( mat, splitExt)
    xy = XY(mat[:, 0], mat[:, 1])
    del mat

    #get peaks from clustering, also get the mean reads per bp to use as global lambda
    rpb = float(xy.number) / (np.max(xy.ys) - np.min(xy.xs))

    if fixybg is not None:
        keybg, mat = parseIxy(fixybg, cut=cut,mcut=mcut)
        if split:
            mat = getSplitMat( mat, splitExt)
        xybg = XY(mat[:, 0], mat[:, 1])
        del mat
        if keybg[0] != key[0]:
            logger.error("%s and %s not the same chromosome. Return." %
                         (fixy, fixybg))
            return
        if len(xybg.xs) == 0:
            logger.error("No PETs found in the control file %s. Return." %
                         fixybg)
            return

    report = "Estimating %s significance with %s, and control %s" % (
        key[0], fixy, fixybg)
    logger.info(report)

    #estimation of peak significance with poisson test
    mpeaks = []
    for peak in tqdm(peaks):
        m = (peak.start + peak.end) / 2
        counts = len(xy.queryPeakBoth(peak.start, peak.end))
        peak.counts = counts
        if peak.counts < max(minPts):
            continue    
        #summit 
        sig = get1DSig(xy, peak.start, peak.end)
        p = np.argmax( sig )
        peak.summit = peak.start + p

        peak.density = counts / 1.0 / peak.length / totalPets * 10.0**9  #RPKM
        peak.poisson_p_value = []
        peak.enrichment_score = []

        us = peak.start - peak.length
        ue = peak.end - peak.length
        ds = peak.start + peak.length
        de = peak.end + peak.length
        ucounts = len(xy.queryPeakBoth(us, ue))
        dcounts = len(xy.queryPeakBoth(ds, de))
        peak.up_down_counts = [ucounts, dcounts]

        #local test
        for ext in exts:
            cs = getPeakNearbyPETs(peak.start, peak.length, xy, ext)
            #median is much sensitive than mean value
            if sen:
                r = min( np.median(cs),np.mean(cs) )
            else:
                r = np.mean(cs)
            if r > 0:
                peak.poisson_p_value.append(
                    max([1e-300, poisson.sf(counts - 1.0, r)]))
                peak.enrichment_score.append(counts / r)
            else:
                peak.poisson_p_value.append(
                    max([1e-300, poisson.sf(counts - 1.0, pseudo)]))
                peak.enrichment_score.append(counts / 1.0)

        #global test
        peak.poisson_p_value.append(
            max([1e-300, poisson.sf(counts - 1.0, rpb * peak.length)]))
        peak.enrichment_score.append(counts / (rpb * peak.length))

        #enrichment score requirement
        if min(peak.enrichment_score) < escut:
            continue

        #get the control counts
        if fixybg is not None and totalControlPets is not None:
            countsbg = len(xybg.queryPeakBoth(peak.start, peak.end))
            peak.control_counts = countsbg
            peak.control_density = countsbg / 1.0 / peak.length / totalControlPets * 10**9

            ucounts = len(xybg.queryPeakBoth(us, ue))
            dcounts = len(xybg.queryPeakBoth(ds, de))
            peak.control_up_down_counts = [ucounts, dcounts]

            if countsbg == 0:
                peak.enrichment_score_vs_control = 100
            else:
                peak.enrichment_score_vs_control = peak.density / peak.control_density
                #peak not significant compared to control, remove
                if peak.enrichment_score_vs_control < ESvsControlCut:
                    continue
            #add nearby local
            peak.control_local_counts = []
            for ext in exts:
                cs = getPeakNearbyPETs(peak.start, peak.length, xybg, ext)
                if sen:
                    r = min( np.median(cs),np.mean(cs) )
                else:
                    r = np.mean(cs)
                peak.control_local_counts.append(r)

        #determine significance
        if max(
                peak.poisson_p_value
        ) <= pcut and peak.length >= lencut:  #if using the min(peak.poisson_p_value), will called more peaks
            peak.significant = 1
        else:
            continue
        mpeaks.append(peak)
    return key, mpeaks


def getFitSf(peaks,cpu=1):
    """
    Get the linear fitting scaling factor using candidate peak nearby regions.
    """
    controlCounts, targetCounts = [], []
    for k, vs in peaks.items():
        for peak in vs:
            for i, c in enumerate(peak.control_up_down_counts):
                controlCounts.append(c / float(peak.length) * 1000)
                targetCounts.append(peak.up_down_counts[i] /
                                    float(peak.length) * 1000)
    sf = estSf(controlCounts, targetCounts, cpu=cpu)
    return sf



def estSigVsControl(key, peaks, sf, pcut=1e-2, pseudo=1, ESvsControlCut=2.0):
    """
    Estimating significance with the control sets.
    @param peaks: list of cLoops2.ds.Peak object
    @param sf: scaling factor
    @param pcut: float, p-value cutoff for poisson test
    @param pseudo: int, pseudo count add to 0 control region, to get reasonable enrichment score and p-value
    @param ESvsControlCut: float, cutoff for enrichment score over control region
    """
    mpeaks = []
    #using poisson test estimate significance vs control
    for peak in peaks:
        cs = [peak.control_counts, pseudo]
        #cs.extend(peak.control_local_counts)
        #normalized backgound counts
        countsbg = max(cs) * sf
        #countsbg = peak.control_counts * sf
        peak.control_scaled_counts = peak.control_counts * sf
        if countsbg > 0:
            peak.poisson_p_value_vs_control = max(
                [1e-300, poisson.sf(peak.counts - 1.0, countsbg)])
        else:
            peak.poisson_p_value_vs_control = 1e-300
        if peak.significant == 1:
            if peak.poisson_p_value_vs_control > pcut:
                continue
            #this will be much strigent for fitting result
            if peak.control_scaled_counts > 0 and peak.counts / peak.control_scaled_counts < ESvsControlCut:
                continue
        mpeaks.append(peak)
    return key, mpeaks



def selSigPeaks(key, peaks):
    """
    Remove overlapped called peaks, keep the more significant one for multiple eps result. 
    """
    report = "Getting the %s most significant peaks for overlapped peaks." % (
        key[0])
    logger.info(report)
    peaks = sortPeaksByLength(peaks)
    #classify peaks to overlapped peaks
    npeaks = []
    skips = set()
    for i in range(len(peaks)):
        if i in skips:
            continue
        n = [peaks[i]]
        for j in range(i + 1, len(peaks)):
            if j in skips:
                continue
            if checkPeakOverlap(peaks[i], peaks[j]):
                skips.add(j)
                n.append(peaks[j])
        npeaks.append(n)
    #final merged peaks
    mpeaks = []
    #sort peak according to significance/density 
    for n in npeaks:
        if len(n) == 1:
            mpeaks.append(n[0])
        else:
            #confirmed through check examples from different eps
            mpeaks = sorted( mpeaks, key=lambda peak: peak.density, reverse=True)
            mpeaks.append(n[0])
    #for a larger peaks overlapped two not-joint peaks, there will be missed. So get these peaks back.
    n = len(mpeaks)
    skips = set()
    while True:
        dels = []
        for i in range(len(peaks)):
            if i in skips:
                continue
            peaka = peaks[i]
            flag = 0
            for peakb in mpeaks:
                if checkPeakOverlap(peaka, peakb):
                    flag = 1
                    skips.add( i )
                    break
            if flag == 0:  #no overlapped significant peaks
                mpeaks.append(peaka)
        if n == len(mpeaks):
            break
        else:
            n = len(mpeaks)
    report = "Calling significant peaks finished for %s, %s non-overlapped peaks. " % (
        key[0], len(mpeaks))
    logger.info(report)
    return mpeaks


def filterPeaksByBonCorr(peaks, pcut=1e-2):
    """
    Bonferroni correction for p-values.
    """
    logger.info("Carrying Bonferroni correction for p-values and further filtering.")
    npeaks = []
    n = len(peaks)
    for peak in tqdm(peaks):
        ps = []
        peak.poisson_p_value = [p * n for p in peak.poisson_p_value]
        ps.extend(peak.poisson_p_value)
        if max(peak.poisson_p_value) > pcut:
            continue
        if hasattr(peak, "poisson_p_value_vs_control"):
            peak.poisson_p_value_vs_control = peak.poisson_p_value_vs_control * n
            if peak.poisson_p_value_vs_control > pcut:
                continue
            ps.append(peak.poisson_p_value_vs_control)
        peak.p_value_mean = len(ps) / sum([1.0 / p for p in ps])
        npeaks.append(peak)
    return npeaks


def callPeaks(
        metaf,
        fout,
        log,
        eps=[100, 200],
        minPts=[3, 5],
        pcut=1e-2,
        cpu=1,
        metabgf=None,
        bgm="ratio",
        pseudo=1,
        sen=False,
        cut=0,
        mcut=1000,
        split=False,
        splitExt=50,
):
    """
    Call peaks main funciton.
    @param metaf: str, petMeta.json file for calling peaks
    @param fout: str, output file prefix
    @param eps: int, eps for blockDBSCAN
    @param minPts: int, minPts for blockDBSCAN
    @param pcut: float, poisson cutoffs for calling significant peaks
    @param cpu: int, cpus to run jobs
    @param metabgf: str, petMeta.json file for control sample
    @param bgm: method to scaling control/input/IgG data, default is normalized with sequencing depth, can be linear fitting
    @param pseudo: int, pseudo count used to control small PETs peaks, 0 is useful for detecting cut&tag sensitive peaks
    @param sen: bool, sensitive mode to call peaks
    @param cut: int, distance >=cut PETs will be used
    @param mcut: int, distance <=mcut PETs will be used
    @param split: bool, whether to split paired-end tags as single-end reads
    @param splitExt: int, extension of the splited single-end read
    """
    global logger
    logger = log
    meta = json.loads(open(metaf).read())
    #get the total PETs
    totalPets = 0
    for key in meta["data"]["cis"].keys():
        f = meta["data"]["cis"][key]["ixy"]
        key, mat = parseIxy(f, cut=cut,mcut=mcut)
        totalPets += mat.shape[0]
    if split:
        totalPets = totalPets * 2

    if sen:
        pseudo = 0

    if metabgf is None:
        #find peaks for files without control
        ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(findSigPeaks)(
            meta["data"]["cis"][key]["ixy"],
            totalPets,
            eps=eps,
            minPts=minPts,
            fixybg=None,
            totalControlPets=None,
            pcut=pcut,
            lencut=max(eps),
            pseudo=pseudo,
            sen=sen,
            cut=cut,
            mcut=mcut,
            split=split,
            splitExt=splitExt,
        ) for key in meta["data"]["cis"].keys())
        peaks = {}
        for d in ds:
            if d is not None:
                peaks[d[0]] = d[1]
    else:
        #using control to furthur determine significance
        metabg = json.loads(open(metabgf).read())
        totalControlPets = 0
        for key in metabg["data"]["cis"].keys():
            f = metabg["data"]["cis"][key]["ixy"]
            key, mat = parseIxy(f, cut=cut,mcut=mcut)
            totalControlPets += mat.shape[0]
        if split:
            totalControlPets = totalControlPets * 2
        #call peaks
        ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(findSigPeaks)(
            meta["data"]["cis"][key]["ixy"],
            totalPets,
            eps=eps,
            minPts=minPts,
            fixybg=metabg["data"]["cis"][key]["ixy"],
            totalControlPets=totalControlPets,
            pcut=pcut,
            lencut=max(eps),
            pseudo=pseudo,
            sen=sen,
            cut=cut,
            mcut=mcut,
            split=split,
            splitExt=splitExt,
        ) for key in meta["data"]["cis"].keys())
        peaks = {}
        for d in ds:
            if d is not None:
                peaks[d[0]] = d[1]

        #estimating scaling factor
        bgsf =  totalPets / float(totalControlPets)
        fitsf = getFitSf( peaks ,cpu=cpu)
        if sen:
            sf = min(bgsf,fitsf)
        else:
            if bgm == "ratio":
                sf = bgsf
            else:
                sf = fitsf
        logger.info( "Estimated scale factor %.3f from linear fitting, meanwhile library ratio is %.3f for target vs control, due to -bgm and -sen, final used scaling factor %.3f" % (fitsf,bgsf, sf))

        #furthur estimate significance vs control
        ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(estSigVsControl)(
            key, peaks[key], sf, pcut=pcut, pseudo=pseudo)
                                  for key in peaks.keys())
        peaks = {}
        for d in ds:
            if d is not None:
                peaks[d[0]] = d[1]

    #for overlapped peaks, furthur get the unique sets
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(selSigPeaks)(
        key,
        peaks[key],
    ) for key in peaks.keys())
    peaks = []
    for d in ds:
        peaks.extend(d)

    #filter by p-values Bonferroni correction
    peaks = filterPeaksByBonCorr(peaks, pcut)

    #write peaks to files
    logger.info("Output %s peaks to %s_peaks.txt and %s_peaks.bed" %
                (len(peaks), fout, fout))
    peaks2txt(peaks, fout + "_peaks.txt")
    peaks2bed(peaks, fout + "_peaks.bed")
