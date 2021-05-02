#!/usr/bin/env python
#--coding:utf-8--
"""
tow.py 
cLoops2 tow algorithm for normalizing ChIP-seq data to reference sample at bp level. 
tmp settings: only remove background noise

Notes
1. log the data, much better, for unfixed size region
2. for bp mean level RPM ,log not needed 
3. should also consider different signal to noise ratio, for peak region
4: two level of normalization: a. normalize to ChIP itself background; b. normalize with scaling 
5: trin GMM model for signals around peaks can have weights of signal and noise, however, the weights are not stable due to the regions selected around peaks, so it is not good choose (signal - wnoise*noise)/wsignal as correction, can be tested

Main function:
1. find background regions 
2. qc for signal/noise ratio and noise level
3. correct the signal based on noise in one sample
4. correct the signal based on noise ratio for treatment sample to reference sample
5. correct the signal regions for treatment sample
"""

#sys library
import time
import sys
import os
import argparse
from glob import glob
from copy import deepcopy
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn.mixture import GaussianMixture as GMM

#cLoops2
from cLoops2.cmat import getBinMean
from cLoops2.utils import getLogger
from cLoops2.settings import *

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Normalize target sample ChIP-seq data signal to reference sample at bp
        level with assumptions of: 1) background noise signal level should be 
        similar and normalized to 0; 2) there are some regions (shared peaks) 
        which share similar level of signal intensities. 

        Input bigWigs file should be normalized as RPM (reads per million) first. 

        The normalization may not suitable for the TF ChIP-seq WT vs KO sample. 

        Example:
        tow.py -br ref.bed -bt tgt.bed -wr ref.bw -wt tgt.bw -o test 
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-br",
        dest="refBed",
        required=True,
        type=str,
        help=
        "The peaks for reference sample in BED format..\n"\
    )
    parser.add_argument(
        "-bt",
        dest="tgtBed",
        required=True,
        type=str,
        help=
        "The peaks for target sample.\n"\
    )
    parser.add_argument("-wr",
                        dest="refBw",
                        required=True,
                        type=str,
                        help="The bigWig file for reference sample.")
    parser.add_argument("-wt",
                        dest="tgtBw",
                        required=True,
                        type=str,
                        help="The bigWig file for target sample.")
    parser.add_argument("-o",
                        dest="fnOut",
                        required=True,
                        type=str,
                        help="Output prefix.")
    parser.add_argument(
        "-ext",
        dest="ext",
        required=False,
        type=int,
        default=10000,
        help=
        "Extension size from peak center to build Gaussian Mixture Model for\n"\
        "classification of background and signal region. Default is 10000bp.\n"\
        "For board peaks such as H3K27me3, should increase the parameter to\n"\
        "include enough nearby regions, such as 50000."
    )
    parser.add_argument("-labelr",
                        dest="refLabel",
                        required=False,
                        type=str,
                        default="ref",
                        help="The label for reference sample. Default is ref.")
    parser.add_argument("-labelt",
                        dest="tgtLabel",
                        required=False,
                        type=str,
                        default="tgt",
                        help="The label for target sample. Default is tgt.")
    op = parser.parse_args()
    return op


def readBed(f):
    rs = []
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        try:
            line[1] = int(line[1])
            line[2] = int(line[2])
        except:
            continue
        rs.append(line[:3])
    return rs


def buildCov(rs):
    """
    Genomic region coverage
    """
    lims = {}
    cov = {}
    for r in rs:
        #coverages
        if r[0] not in cov:
            cov[r[0]] = set()
        cov[r[0]].update(range(r[1], r[2] + 1))
        #range limitations
        if r[0] not in lims:
            lims[r[0]] = [r[1], r[2]]
        if r[1] < lims[r[0]][0]:
            lims[r[0]][0] = r[1]
        if r[2] > lims[r[0]][1]:
            lims[r[0]][1] = r[2]
    return cov, lims


def getRegion(cov, margin=1):
    """
    Get regions from coverage
    """
    rs = []
    for c, s in cov.items():
        s = list(s)
        s.sort()
        i = 0
        while i < len(s) - 1:
            for j in range(i + 1, len(s)):
                if s[j] - s[j - 1] > margin:
                    break
                else:
                    continue
            start = s[i]
            end = s[j - 1]
            i = j  #update search start
            rs.append([c, start, end])
    return rs


def checkBgOverlaps(c, s, e, cov, lims):
    if s < lims[c][0] or e > lims[c][1]:
        return True
    for i in range(s, e + 1):
        if i in cov[c]:
            return True
    return False


def getFgBgs(refPeaks, tgtPeaks, ext=10):
    """
    Get the background regions.
    #coverage, require no overlaps for background regions
    """
    refCov, refLims = buildCov(refPeaks)
    tgtCov, tgtLims = buildCov(tgtPeaks)
    coCov, allCov, lims = {}, {}, {}
    for c in refCov.keys():
        minv = min(refLims[c][0], tgtLims[c][0])
        maxv = max(refLims[c][1], tgtLims[c][1])
        lims[c] = [minv, maxv]
        a = deepcopy(refCov[c])
        a.update(deepcopy(tgtCov[c]))
        allCov[c] = a
        s = refCov[c].intersection(tgtCov[c])
        coCov[c] = s
    fgs = getRegion(coCov)  #shared peaks regions
    alls = getRegion(allCov)  #all peak regions
    bgs = []
    for r in alls:
        d = r[2] - r[1]
        c = r[0]
        s = r[1] - ext * d
        e = r[2] - ext * d
        if checkBgOverlaps(c, s, e, allCov, lims):
            continue
        bgs.append([c, s, e])
        s = r[1] + ext * d
        e = r[2] + ext * d
        if checkBgOverlaps(c, s, e, allCov, lims):
            continue
        bgs.append([c, s, e])
    return fgs, bgs


def getBwSig(rs, f, bins=100):
    """
    Get region signal from bigWig file. 
    """
    bw = pyBigWig.open(f)
    s = []
    for r in rs:
        try:
            ns = bw.values(r[0], r[1], r[2])
        except:
            continue
        ns = np.nan_to_num(ns)
        if len(ns) < bins:
            continue
        ns = getBinMean(ns, bins)
        s.append(ns)
    return np.array(s)


def getQc(fgs,
          bgs,
          refBw,
          tgtBw,
          fnOut,
          refLabel="ref",
          tgtLabel="tgt",
          bins=100):
    """
    Quality control: 
    1. noise compare
    2. signal to noise ratio compare
    """
    fgRef = getBwSig(fgs, refBw, bins=bins)
    fgRef = fgRef.mean(axis=0)
    fgTgt = getBwSig(fgs, tgtBw, bins=bins)
    fgTgt = fgTgt.mean(axis=0)
    bgRef = getBwSig(bgs, refBw, bins=bins)
    bgRef = bgRef.mean(axis=0)
    bgTgt = getBwSig(bgs, tgtBw, bins=bins)
    bgTgt = bgTgt.mean(axis=0)
    bgfc = bgRef.mean() / bgTgt.mean()
    fgfc = fgRef.mean() / fgTgt.mean()
    fig, axs = pylab.subplots(1, 3, figsize=(6, 2), sharex=True)
    axs = axs.reshape(-1)
    x = np.arange(bins)
    ax = axs[0]
    ax.plot(x, bgRef, label=refLabel)
    ax.plot(x, bgTgt, label=tgtLabel)
    ax.set_xlabel("bins")
    ax.set_ylabel("ChIP-seq mean signals, RPM")
    ax.set_title("background region\nsf(%s->%s):%.3f" %
                 (tgtLabel, refLabel, bgfc))
    ax.legend()
    #signal level
    ax = axs[1]
    ax.plot(x, fgRef, label=refLabel)
    ax.plot(x, fgTgt, label=tgtLabel)
    ax.set_ylabel("ChIP-seq mean signals, RPM")
    ax.set_xlabel("bins")
    ax.set_title("peak region\nsf(%s->%s):%.3f" % (tgtLabel, refLabel, fgfc))
    ax.legend()
    #signal to noise ratio
    refSN = fgRef / bgRef
    tgtSN = fgTgt / bgTgt
    ax = axs[2]
    ax.plot(x, refSN, label=refLabel)
    ax.plot(x, tgtSN, label=tgtLabel)
    ax.set_xlabel("bins")
    ax.set_ylabel("Signal to noise ratio")
    ax.legend()
    ax.set_title("Signal to noise")
    pylab.tight_layout()
    pylab.savefig(fnOut + "_1_qc.pdf")
    return fgRef, fgTgt, bgRef, bgTgt


def estFit(bgs,
           refBw,
           tgtBw,
           bgRef,
           bgTgt,
           sf,
           fnOut,
           refLabel="ref",
           tgtLabel="tgt",
           bins=2):
    """
    Estimate linear fiting between samples.
    """
    refS = getBwSig(bgs, refBw, bins=bins)
    refS = pd.Series(refS.reshape(-1))
    tgtS = getBwSig(bgs, tgtBw, bins=bins)
    tgtS = pd.Series(tgtS.reshape(-1))

    fig, axs = pylab.subplots(1, 2, figsize=(5, 2))
    axs = axs.reshape(-1)

    #signal conversion
    refS = refS - bgRef.mean()
    tgtS = (tgtS - bgTgt.mean()) * sf
    s = refS[refS > 0].index
    s = tgtS[s]
    s = s[s > 0].index
    refS = refS[s]
    tgtS = tgtS[s]
    #log transformation
    refS = np.log2(refS)
    tgtS = np.log2(tgtS)
    ax = axs[0]
    sns.kdeplot(refS, label=refLabel, ax=ax, fill=True)
    sns.kdeplot(tgtS, label=tgtLabel, ax=ax, fill=True)
    ax.legend()
    ax.set_xlabel("log2(RPM)")
    ax.set_title("signal distribution")

    #distribution match
    ax = axs[1]
    #tgtSc = (tgtS - tgtS.mean())/tgtS.std()*refS.std() + refS.mean()
    alpha = refS.std() / tgtS.std()
    beta = refS.mean() - alpha * tgtS.mean()
    tgtSc = tgtS * alpha + beta

    m = refS - tgtSc
    a = (refS + tgtSc) / 2
    ax.scatter(a, m, s=0.1)
    if beta > 0:
        ax.set_title(
            "after correction M~A PCC:%.3f\nlog2(%s)=%.3flog2(%s)+%.3f" %
            (m.corr(a), refLabel, alpha, tgtLabel, beta))
    else:
        ax.set_title(
            "after correction M~A PCC:%.3f\nlog2(%s)=%.3flog2(%s)%.3f" %
            (m.corr(a), refLabel, alpha, tgtLabel, beta))
    ax.set_xlabel("A, (log2(%s)+log2(%s))/2)" % (refLabel, tgtLabel))
    ax.set_ylabel("M, log2(%s)-log2(%s)" % (refLabel, tgtLabel))

    pylab.tight_layout()
    pylab.savefig(fnOut + "_2_fgSignalConversion.pdf")
    return [alpha, beta]


def corrSig(ss, noise=None, sf=None, sf2=None, trim=False):
    """
    ss: signal matrix
    noise: value, random noise
    sn: signal noise ratio
    sf: scaling fitting, if none, do not scaling
    """
    ns = []
    ss = pd.DataFrame(ss)
    for t in ss.itertuples():
        t = np.array(t[1:])
        if noise is not None:
            t = t - noise
        if trim:
            t[t < 0] = 0
        if sf is not None:
            t = t * sf
        if sf2 is not None:
            t = 2**(np.log2(t) * sf2[0] + sf2[1])
        ns.append(t)
    return pd.DataFrame(ns)


def checkCorrSig(rs,
                 bgs,
                 refBw,
                 tgtBw,
                 bgRef,
                 bgTgt,
                 sf,
                 sf2,
                 fnOut,
                 refLabel="ref",
                 tgtLabel="tgt",
                 bins=100):
    """
    Correct signal.
    """
    fgRef = corrSig(getBwSig(rs, refBw, bins=bins), bgRef,
                    trim=True).mean(axis=0)
    bgRef = corrSig(getBwSig(bgs, refBw, bins=bins), bgRef,
                    trim=False).mean(axis=0)

    fgTgt = corrSig(getBwSig(rs, tgtBw, bins=bins), bgTgt, sf, sf2,
                    trim=True).mean(axis=0)
    bgTgt = corrSig(getBwSig(bgs, tgtBw, bins=bins), bgTgt, sf,
                    trim=False).mean(axis=0)

    fig, axs = pylab.subplots(1, 2, figsize=(6, 3), sharex=True)
    axs = axs.reshape(-1)
    #noise level
    x = np.arange(bins)
    ax = axs[0]
    ax.plot(x, bgRef, label=refLabel)
    ax.plot(x, bgTgt, label=tgtLabel)
    ax.set_ylim([-0.001, 0.001])
    ax.set_xlabel("bins")
    ax.set_ylabel("ChIP-seq mean signals, RPM")
    ax.set_title("background region")
    ax.legend()
    #signal level
    ax = axs[1]
    ax.plot(x, fgRef, label=refLabel)
    ax.plot(x, fgTgt, label=tgtLabel)
    ax.set_ylabel("ChIP-seq mean signals, RPM")
    ax.set_xlabel("bins")
    ax.set_title("peak region")
    ax.legend()
    pylab.tight_layout()
    pylab.savefig(fnOut + "_3_CorrectSignalQc.pdf")


def getGmm(fg, bg, bw, ax, title, ext=5000, bins=5):
    """
    Train one GMM. Data should log2 first
    """
    rs = []
    for r in fg:
        center = int((r[1] + r[2]) / 2)
        nr = [r[0], center - ext, center + ext]
        rs.append(nr)
    #foreground signal
    fgS = getBwSig(fg, bw, bins=bins)
    fgS = fgS.reshape(-1)
    fgS = np.log2(fgS[fgS > 0])
    #background signal
    bgS = getBwSig(bg, bw, bins=bins)
    bgS = bgS.reshape(-1)
    bgS = np.log2(bgS[bgS > 0])
    #peak nearby region
    mixS = getBwSig(rs, bw, bins=bins)
    mixS = mixS.reshape(-1)
    mixS = np.log2(mixS[mixS > 0])
    #plot
    sns.kdeplot(bgS, label="background", fill=False, ax=ax)
    sns.kdeplot(fgS, label="peaks", fill=False, ax=ax)
    sns.kdeplot(mixS, label="peak and nearby", fill=False, ax=ax)
    #train gmm
    gmm = GMM(n_components=2,
              covariance_type="full",
              random_state=123,
              means_init=[[np.mean(fgS)],
                          [np.mean(bgS)]]).fit([[v] for v in mixS])
    ms = gmm.means_.reshape(-1)
    ws = gmm.weights_
    ax.legend()
    ax.set_xlabel("log2(RPM)")
    ax.set_title(title + "\n" + "means:%.3f, %.3f\nweights:%.3f, %.3f" %
                 (ms[0], ms[1], ws[0], ws[1]))
    #correct gmm predict targets, 0 as noise , 1 as signal
    if ms[0] > ms[1]:
        cs = {0: 1, 1: 0}
    else:
        cs = {0: 0, 1: 1}
    return gmm, cs


def trainGmm(refPeaks,
             tgtPeaks,
             bgs,
             refBw,
             tgtBw,
             fnOut,
             refLabel="ref",
             tgtLabel="tgt",
             bins=5):
    """
    Train GMM to classify background regions or signal regions
    """
    fig, axs = pylab.subplots(1, 2, figsize=(6, 3), sharex=True, sharey=True)
    axs = axs.reshape(-1)

    refGmm, refCs = getGmm(refPeaks, bgs, refBw, axs[0], refLabel)
    tgtGmm, tgtCs = getGmm(tgtPeaks, bgs, tgtBw, axs[1], tgtLabel)

    pylab.tight_layout()
    pylab.savefig(fnOut + "_4_GMM.pdf")
    return refGmm, refCs, tgtGmm, tgtCs


def normRefBw(bw, gmm, gmmCs, noise, fnOut):
    """
    Normalize reference bigWig files. Only remove background noise. 
    """
    bwi = pyBigWig.open(bw)
    with open(fnOut + ".bdg", "w") as fo:
        for chrom, size in bwi.chroms().items():
            #get the singal for whole chromosome
            ss = bwi.intervals(chrom)
            for s in ss:
                v = s[-1]
                if v == 0:
                    continue
                #gmm was trained with log2 data
                t = gmmCs[gmm.predict([[np.log2(v)]])[0]]
                if t == 0:
                    continue
                v = v - noise
                if v < 0:
                    continue
                line = [chrom, s[0], s[1], "%.5f" % v]
                fo.write("\t".join(list(map(str, line))) + "\n")
    bwi.close()


def normTgtBw(bw, gmm, gmmCs, noise, sf, sf2, fnOut):
    """
    Normalize target sample bigWig file.
    """
    bwi = pyBigWig.open(bw)
    with open(fnOut + ".bdg", "w") as fo:
        for chrom, size in bwi.chroms().items():
            #add header
            ss = bwi.intervals(chrom)
            for s in ss:
                v = s[-1]
                if v == 0:
                    continue
                #gmm was trained with log2 data
                t = gmmCs[gmm.predict([[np.log2(v)]])[0]]
                if t == 0:
                    continue
                v = v - noise
                if v < 0:
                    continue
                #first scaling factor, obtained with not log2 data
                v = v * sf
                if t == 1:  #signal, further normalization
                    v = 2**(np.log2(v) * sf2[0] + sf2[1])
                line = [chrom, s[0], s[1], "%.5f" % v]
                fo.write("\t".join(list(map(str, line))) + "\n")
    bwi.close()


def main():
    op = help()

    #step 0 parameters check
    report = "python tow.py -br {refBed} -bt {tgtBed} -wr {refBw} -wt {tgtBw} -o {fnOut} -ext {ext} -labelr {refLabel} -labelt {tgtLabel}".format(
        refBed=op.refBed,
        tgtBed=op.tgtBed,
        refBw=op.refBw,
        tgtBw=op.tgtBw,
        fnOut=op.fnOut,
        ext=op.ext,
        refLabel=op.refLabel,
        tgtLabel=op.tgtLabel)
    logger.info(report)

    #step 1 read bed
    refPeaks = readBed(op.refBed)
    tgtPeaks = readBed(op.tgtBed)

    #step 2 generate background regions
    fgs, bgs = getFgBgs(refPeaks, tgtPeaks)
    logger.info(
        "[%s] Step0: ref sample peaks: %s; tgt sample peaks: %s; shared: %s; background regions: %s"
        % (op.fnOut, len(refPeaks), len(tgtPeaks), len(fgs), len(bgs)))

    #step 3 qc for signal to noise ratio and noise level
    logger.info(
        "[%s] Step1: Initial QC for background noise level and signal-to-noise ratio."
        % (op.fnOut))
    fgRef, fgTgt, bgRef, bgTgt = getQc(fgs,
                                       bgs,
                                       op.refBw,
                                       op.tgtBw,
                                       op.fnOut,
                                       refLabel=op.refLabel,
                                       tgtLabel=op.tgtLabel)
    #scaling factor for background region
    sf = bgRef.mean() / bgTgt.mean()

    #step 4 estimate sample-wise fitting
    #scaling factor for singla region
    logger.info("[%s] Step2: Estimating scaling factor for signal regions." %
                (op.fnOut))
    sf2 = estFit(fgs, op.refBw, op.tgtBw, bgRef, bgTgt, sf, op.fnOut,
                 op.refLabel, op.tgtLabel)

    #step 5 check correction signal
    logger.info("[%s] Step3: Checking corrected signal." % (op.fnOut))
    checkCorrSig(fgs, bgs, op.refBw, op.tgtBw, bgRef, bgTgt, sf, sf2, op.fnOut,
                 op.refLabel, op.tgtLabel)

    #step 6 train GMM with fg and bg data for classifiy fg and bg regions
    logger.info(
        "[%s] Step4: Building Gaussian Mixture Model for classfication of background and signal regions"
        % (op.fnOut))
    refGmm, refCs, tgtGmm, tgtCs = trainGmm(refPeaks,tgtPeaks, bgs, op.refBw, op.tgtBw,
                                            op.fnOut, op.refLabel, op.tgtLabel)

    #step 7 performe normalization to bigWig files
    noiseRef = bgRef.mean()
    noiseTgt = bgTgt.mean()
    logger.info(
        "[%s] Step5: Normaling reference sample with background noise." %
        (op.fnOut))
    normRefBw(op.refBw, refGmm, refCs, noiseRef, op.fnOut + "_" + op.refLabel)
    logger.info("[%s] Step6: Normaling target sample." % (op.fnOut))
    normTgtBw(op.tgtBw, tgtGmm, tgtCs, noiseTgt, sf, sf2,
              op.fnOut + "_" + op.tgtLabel)

    logger.info("[%s] analysis finished." % (op.fnOut))


if __name__ == "__main__":
    start_time = datetime.now()
    main()
    usedtime = datetime.now() - start_time
    sys.stderr.write("Process finished. Used CPU time: %s Bye!\n" % usedtime)
