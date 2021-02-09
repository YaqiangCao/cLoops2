#!/usr/bin/env python3
#--coding:utf-8 --
"""
1. Insulation score does not work well for Hi-Trac/Trac-looping data.
2. z-score normalization is much better than log2(s/mean), much stable
3. 20k similar to 10k and 5k, 1k will not improve the quality of domains called and score showed, much time consuming.
4. fix such as 10k and 5k, fine tune naerby 500k to 250k, will affect a lot. Better just use 10k and 500k. 
5. Variants of insulation score, may not work as stable as SS.
6. when caculate SS, not remove <0 as 0, will cause strange results.
7. if multiple parameters given, call nest domains?
8. obs/exp matrix than correlation, not work, by all means

2020-03-08: update density
2020-09-13: update multiple window size added
2020-09-22: try to improve efficiency of function calcSS, by getting the whole chromosome contact matrix as sparse matrix. Seems sparse matrix will auto occupy multipe CPUs.
2020-11-18: going to integrate insulation score for hic
"""

#sys
import json
from copy import deepcopy

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#cLoops2
from cLoops2.ds import XY,Domain
from cLoops2.io import parseIxy, doms2txt, doms2bed
from cLoops2.cmat import getObsMat, xy2dict, dict2mat
from cLoops2.settings import *




def calcSS(f, bs=20000, winSize=500000, cut=0,mcut=-1,hic=False):
    """
    Calculation of correlation matrix insulation score, output as .bedGraph file.
    @param bs: bin size
    @param winSize: sliding matrix width half size
    @param cut: distance cutoff for PETs
    """
    key, mat = parseIxy(f, cut=cut,mcut=mcut)
    matstart = np.min(mat)
    matend = np.max(mat)
    start = matstart + winSize
    end = matend - winSize
    bins = int((end - start) / bs)
    #convert to sparse contact matrix
    mat = xy2dict(mat, s=matstart, e=matend, r=bs)
    mat = dict2mat(mat)
    print(
        "caculating from %s to %s of %s bins for segregation score with bin size of %s and window size of %s"
        % (start, end, bins, bs,winSize))
    rs = []
    ss = []
    for i in tqdm(range(bins)):
        x = start + i * bs
        s = x - winSize
        e = x + winSize
        #releative position in contact matrix
        s = int( (s - matstart)/bs )
        e = int( (e - matstart)/bs ) +1
        nmat = mat[s:e,s:e].toarray()
        #previous 
        #nmat = getObsMat(mat, s, e, bs)
        nmat = np.log2(nmat + 1)
        nmat = np.corrcoef(nmat)
        nmat = np.nan_to_num(nmat)
        nmat = nmat[int(nmat.shape[0] / 2) + 1:, :int(nmat.shape[1] / 2)]
        if hic == False:
            nmat[nmat < 0] = 0
        s = nmat.mean()
        ss.append(s)
        r = [key[0], x, x + bs]
        rs.append(r)
    ss = np.array(ss)
    ss = (ss - np.mean(ss)) / np.std(ss)
    for i, r in enumerate(rs):
        r.append(ss[i])
    return key, rs


def writeSS2Bdg(rs, fout):
    """
    Write segragation score as bedGraph file.
    """
    with open(fout, "w") as fo:
        for r in rs:
            fo.write("\t".join(list(map(str, r))) + "\n")


def callDom(key, rs, bs, ws, cut=0, mcut=-1, lencut=10):
    """
    Call domain based on caculated segragation score.
    """
    #find domains
    doms = []
    i = 0
    while i < len(rs):
        if rs[i][-1] > cut:
            j = i + 1
            p = -1
            while j < len(rs):
                if rs[j][-1] > cut:
                    p = j
                    j += 1
                else:
                    break
            if p > i and rs[p][2] - rs[i][1] > lencut * bs:
                dom = Domain()
                dom.chrom = rs[i][0]
                dom.start = rs[i][1]
                dom.end = rs[p][2]
                dom.length = dom.end - dom.start
                dom.ss = np.mean([t[-1] for t in rs[i:p + 1]])
                dom.bs = bs
                dom.ws = ws
                doms.append(dom)
            i = j
        else:
            i = i + 1
    return "-".join(key), doms


def compDoms(doma,domb,lrcut=0.9):
    """
    Compare if is quite close same domains.
    If quite close, whether use doma to replace domb.
    """
    if doma.chrom != domb.chrom:
        return False, None
    #overlapped domains
    if domb.start <= doma.start <= domb.end or domb.start <= doma.end <= domb.end or doma.start <= domb.start <= doma.end or doma.start <= domb.end <= doma.end:
        start = max(doma.start,domb.start)
        end = min(doma.end,domb.end)
        length = max(doma.length,domb.length)
        if (end-start)/length > lrcut :
            if doma.ss > domb.ss:
                return True, doma
            return True, None
    return False, None
 

def combineDoms(doms, doms2, lrcut=0.9):
    """
    Combine domains.
    """
    #doms binsize is smaller than doms2
    for key in doms2.keys():
        if key not in doms:
            doms[key] = doms2[key]
        else:
            #add no overlapped
            for doma in doms2[key]:
                flag = False
                for i, domb in enumerate(doms[key]):
                    flag2,n = compDoms(doma,domb,lrcut) #if highly overlapped and similar, do not add the new record
                    if flag2:
                        flag = True
                        if n is not None:
                            doms[key][i] = doma #replace 
                        break
                    else:
                        continue
                #no overlapped or almost the same
                if flag == False:
                    doms[key].append(doma)
    return doms


def quantifyDom(f, doms, tot,cut=0,mcut=-1,hic=False):
    """
    Quantify domains
    """
    key, mat = parseIxy(f, cut=cut,mcut=mcut)
    if mat.shape[0] == 0:
        print(
            "No PETs found in %s maybe due to distance cutoff for PET > %s <%s." %
            (fixy, cut,mcut))
        return None, None
    xy = XY(mat[:, 0], mat[:, 1])
    print("Quantify %s domains from %s" % (len(doms), key))
    ndoms = []
    for dom in tqdm(doms):
        t = xy.queryPeak(dom.start, dom.end)
        b = xy.queryPeakBoth(dom.start, dom.end)
        n = t.difference(b)
        if len(n) > 0:
            e = len(b) / float(len(n))
        else:
            e = 100
        if hic==False and  e < 1:
            continue 
        dom.totalPETs = len(b) + len(n)
        dom.withinDomainPETs = len(b)
        dom.enrichmentScore = e
        dom.density = dom.withinDomainPETs/float(tot)/float(dom.length)*10**9
        ndoms.append( dom )
    return key, ndoms



def callDomains(
        metaf,
        fout,
        logger,
        bs=[5000, 10000],
        ws=[500000],
        cut=0,
        mcut=-1,
        cpu=1,
        hic=False
):
    """
    Call domains main funciton.
    @param metaf: str, petMeta.json file for calling peaks
    @param fout: str, output file prefix
    @param bs: list of int, bin size for calling domains
    @param ws: list of int, window size for caculating segregation score
    """
    meta = json.loads(open(metaf).read())
    doms = {}  #candidate doamins
    bs.sort()
    tot = meta["Unique PETs"]
    for binSize in bs:
        for winSize in ws:
            #caculating scores
            ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(calcSS)(
                meta["data"]["cis"][key]["ixy"],
                bs=binSize,
                winSize=winSize,
                cut=cut,
                mcut=mcut,
                hic=hic,
            ) for key in meta["data"]["cis"].keys())
            rs = []
            for d in ds:
                rs.extend(d[1])
            writeSS2Bdg(rs, fout + "_domains_SS_binSize%sk_winSize%sk.bdg" % (binSize / 1000,winSize/1000))
            #call domains
            ds = {d[0]: d[1] for d in ds}
            ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(callDom)(key, ds[key], binSize,winSize)
                                      for key in ds.keys())
            doms_2 = {d[0]: d[1] for d in ds}
            doms = combineDoms(doms, doms_2)
    #furthur quantify domains
    ds = Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(quantifyDom)(
        meta["data"]["cis"][key]["ixy"],
        doms[key],
        tot,
        cut,
        mcut,
        hic
    ) for key in doms.keys())
    doms = []
    for d in ds:
        doms.extend(d[1])
    #output domains
    doms2txt(doms, fout + "_domains.txt")
    doms2bed(doms, fout + "_domains.bed")
