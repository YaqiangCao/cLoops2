#!/usr/bin/env python
#--coding:utf-8--
"""
Quantify peaks, loops and domains.
2020-03-09: change quantify loop density with total library size
2020-03-12: update 0 counts peaks/loops significant to totally not significant
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
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from scipy.stats import hypergeom, binom, poisson

#cLoops2
from cLoops2.ds import XY, Loop, Peak
from cLoops2.io import parseTxt2Loops, parseTxt2Domains, ixy2pet, doms2txt, parseIxy,parseBed2Peaks
from cLoops2.settings import *
from cLoops2.callPeaks import getPeakNearbyPETs
from cLoops2.callDomains import calcSS,writeSS2Bdg
from cLoops2.callCisLoops import getPerRegions, estAnchorSig

### peaks quantification releated funcitons
def _quantPeaks(key,peaks,fixy,tot,cut=0, mcut=-1,pext=0,exts=[5,10]):
    key2, mat = parseIxy(fixy, cut=cut,mcut=mcut)
    xy = XY(mat[:, 0], mat[:, 1])
    rpb = float(xy.number) / (np.max(xy.ys) - np.min(xy.xs))
    print("%s \t quantify %s candidate peaks in %s." %
          (datetime.now(), len(peaks), key))
    for peak in tqdm(peaks):
        counts = len(xy.queryPeakBoth(peak.start-pext, peak.end+pext))
        peak.counts = counts
        if peak.counts == 0:
            peak.density = 0
            peak.poisson_p_value = [1,1,1]
            peak.enrichment_score = [0,0,0]
        else:
            peak.density = counts / 1.0 / peak.length / tot * 10.0**9  #RPKM
            peak.poisson_p_value = []
            peak.enrichment_score = []
            #local test
            for ext in exts:
                cs = getPeakNearbyPETs(peak.start, peak.length, xy, ext)
                r = np.median(cs)
                if r > 0:
                    peak.poisson_p_value.append( max([1e-300, poisson.sf(counts - 1.0, r)]))
                    peak.enrichment_score.append(counts / r)
                else:
                    peak.enrichment_score.append(100)
                    peak.poisson_p_value.append(1e-300)
            #global test
            peak.poisson_p_value.append( max([1e-300, poisson.sf(counts - 1.0, rpb * peak.length)]))
            peak.enrichment_score.append(counts / (rpb * peak.length))
    return key,peaks


def _peaks2txt(peaks, fout):
    """
    Converting list of cLoops2.ds.Peaks objects into txt file.
    """
    with open(fout, "w") as fo:
        header = [
            "peakId", "chr", "start", "end", "length", "counts", "RPKM",
            "enrichmentScore", "poissonPvalue"
        ]
        fo.write("\t".join(header) + "\n")
        for peak in peaks:
            line = [
                peak.id, peak.chrom, peak.start, peak.end, peak.length,
                peak.counts, peak.density, peak.enrichment_score,
                peak.poisson_p_value
            ]
            fo.write("\t".join(list(map(str, line))) + "\n")


def quantPeaks(
        predir,
        peakf,
        output,
        logger,
        cut=0,
        mcut=-1,
        cpu=1,
):
    """
    Quantification of peaks.
    """
    #read in peaks orgainzed by chromosomes
    peaks = parseBed2Peaks(peakf)
    #meta data
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(peaks.keys())))
    tot = meta["Unique PETs"]
    #get the data
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(_quantPeaks)(
        key,
        peaks[key],
        meta["data"]["cis"][key]["ixy"],
        tot,
        cut=cut,
        mcut=mcut,
    ) for key in keys)
    peaks = []
    for d in ds:
        peaks.extend(d[1])
    _peaks2txt(peaks, output + "_peaks.txt")
 


### loops quantification releated functions
def _quantLoops(key, loops, fixy, tot, pcut=0, mcut=-1, pseudo=1):
    """
    Estimate the loop density and statstical significance for one chromosomal.
    @param key: str, such as chr21-chr21
    @param loops: list of Loop object
    @param fixy: cLoops2 pre generated .ixy file
    """
    tot = float(tot)
    xy = ixy2pet(fixy, cut=pcut,mcut=mcut)
    N = xy.number
    print("%s \t quantify %s candidate loops in %s." %
          (datetime.now(), len(loops), key))
    nloops = []
    for loop in tqdm(loops):
        ra, rb, rab = xy.queryLoop(loop.x_start, loop.x_end, loop.y_start,
                                   loop.y_end)
        ra, rb, rab = len(ra), len(rb), len(rab)
        #make sure the anchor are significant
        px, esx = estAnchorSig(xy, loop.x_start, loop.x_end)
        py, esy = estAnchorSig(xy, loop.y_start, loop.y_end)
        loop.x_peak_poisson_p_value = px
        loop.x_peak_es = esx
        loop.y_peak_poisson_p_value = py
        loop.y_peak_es = esy
        lowerra, lowerrb, lowerrab = xy.queryLoop(
            loop.x_start - (loop.x_end - loop.x_start), loop.x_start,
            loop.y_start - (loop.y_end - loop.y_start), loop.y_start)  #p2ll
        loop.P2LL = float(rab) / max(len(lowerrab), pseudo)
        loop.ra = ra
        loop.rb = rb
        loop.rab = rab
        if rab > 0:
            #hypergeometric p-value
            hyp = max([1e-300, hypergeom.sf(rab - 1.0, N, ra, rb)])
            #print(ra,rb,rab, np.log10(loop.distance),N,n)
            #start caculate the permutated background
            nas, nbs = getPerRegions(loop, xy)
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
                    else:
                        rabs.append(0)
                        nbps.append(0.0)
            rabs, nbps = np.array(rabs), np.array(nbps)
            if np.median(rabs) > 0:
                mrabs = float(np.median(rabs))
            else:
                mrabs = pseudo
            if np.median(nbps) > 0:
                mbps = np.median(nbps)
            else:
                mbps = 1e-10
            #print(mrabs,mbps,loop.rab,loop.x_end-loop.x_start, loop.y_end-loop.y_start, N)
            #local fdr
            if len(rabs) > 0:
                fdr = len(rabs[rabs > rab]) / float(len(rabs))
            else:
                fdr = 0.0
            #enrichment score
            es = rab / mrabs
            #simple possion test
            pop = max([1e-300, poisson.sf(rab - 1.0, mrabs)])
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
        else:
            loop.FDR = 1
            loop.ES = 0
            loop.density = 0
            loop.hypergeometric_p_value = 1
            loop.poisson_p_value = 1
            loop.binomial_p_value = 1
        nloops.append(loop)
    return key, nloops


def _loops2txt(loops, fout):
    """
    Converting list of cLoops2.ds.loops objects into txt file.
    """
    with open(fout, "w") as fo:
        header = [
            "loopId", "chrA", "startA", "endA", "chrB", "startB", "endB",
            "distance(bp)", "centerA", "centerB", "readsA", "readsB", "cis",
            "PETs", "density", "enrichmentScore", "P2LL", "FDR",
            "binomialPvalue", "hypergeometricPvalue", "poissonPvalue",
            "poissonPvaluePeakA", "poissonPvaluePeakB"
        ]
        fo.write("\t".join(header) + "\n")
        for i, loop in enumerate(loops):
            line = [
                loop.id,
                loop.chromX,
                loop.x_start,
                loop.x_end,
                loop.chromY,
                loop.y_start,
                loop.y_end,
                loop.distance,
                loop.x_center,
                loop.y_center,
                loop.ra,
                loop.rb,
                loop.cis,
                loop.rab,
                loop.density,
                loop.ES,
                loop.P2LL,
                loop.FDR,
                loop.binomial_p_value,
                loop.hypergeometric_p_value,
                loop.poisson_p_value,
                loop.x_peak_poisson_p_value,
                loop.y_peak_poisson_p_value,
            ]
            fo.write("\t".join(list(map(str, line))) + "\n")


def quantLoops(
        predir,
        loopf,
        output,
        logger,
        cut=0,
        mcut=-1,
        cpu=1,
):
    """
    Quantification of loops.
    """
    loops = parseTxt2Loops(loopf, cut=0)
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    tot = meta["Unique PETs"]
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(loops.keys())))
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(_quantLoops)(
        key,
        loops[key],
        meta["data"]["cis"][key]["ixy"],
        tot,
        pcut=cut,
        mcut=mcut,
    ) for key in keys)
    loops = []
    for d in ds:
        loops.extend(d[1])
    _loops2txt(loops, output + "_loops.txt")


### domains quantification releated functions
def _quantDomains(key, domains, fixy,tot, bs=10000, ws=500000, cut=0,mcut=-1):
    """
    Quantify domains
    """
    key, mat = parseIxy(fixy, cut=cut,mcut=mcut)
    if mat.shape[0] == 0:
        print(
            "No PETs found in %s maybe due to distance cutoff for PET > %s <%s." %
            (fixy, cut,mcut))
        return None, None
    xy = XY(mat[:, 0], mat[:, 1])
    keys, rs = calcSS(fixy, bs=bs, winSize=ws, cut=cut)
    print("Quantify %s domains from %s" % (len(domains), key))
    for dom in tqdm(domains):
        t = xy.queryPeak(dom.start, dom.end)
        b = xy.queryPeakBoth(dom.start, dom.end)
        n = t.difference(b)
        if len(n) > 0:
            e = len(b) / float(len(n))
        else:
            e = np.inf
        dom.totalPETs = len(b) + len(n)
        dom.withinDomainPETs = len(b)
        dom.enrichmentScore = e
        dom.density = dom.withinDomainPETs / float(tot) / float(
            dom.length) * 10**9
        dom.bs = bs
        dom.ws = ws
        ps = int((dom.start - np.min(mat) - ws) / bs)
        if ps < 0:
            ps = 0
        pe = int((dom.end - np.min(mat) - ws) / bs) #exact the same result to callDomains
        if pe <0:
            pe = 0
        dom.ss = np.mean([t[-1] for t in rs[ps:pe]])
    return key, domains, rs


def quantDomains(
        predir,
        domainf,
        output,
        logger,
        bs=10000,
        ws=500000,
        cut=0,
        mcut=-1,
        cpu=1,
        bdg=False,
):
    """
    Quantification of domains.
    """
    #domains
    domains = parseTxt2Domains(domainf)
    #pre data
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys = list(set(keys).intersection(set(domains.keys())))
    tot = meta["Unique PETs"]

    #get
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(_quantDomains)(
        key, domains[key], meta["data"]["cis"][key]["ixy"],tot, bs, ws, cut,mcut)
                              for key in keys)
    domains = []
    for d in ds:
        domains.extend(d[1])
    doms2txt(domains, output + "_domains.txt")
    
    #output the score as bedGraph
    if bdg:
        rs = []
        for d in ds:
            rs.extend(d[-1])
        writeSS2Bdg(rs, output + "_SS.bdg")


