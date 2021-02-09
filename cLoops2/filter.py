#!/usr/bin/env python
#--coding:utf-8--
"""
filter.py
cLoops2 PETs filtering related code.

2020-01-29: --invert-match added.
2020-02-17: filter singleton PETs added
2020-02-19: sample PETs added
2020-04-06: filterPETsbyLoops updated as requiring both ends in loops
2020-08-24: update np.random.choice, its default replace parameter is True, there fore the sampling not well for tot<true unique reads.
"""

__author__ = "CAO Yaqiang"
__date__ = "2020-01-28"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import json
from glob import glob
from datetime import datetime

#3rd
import joblib
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed

#cLoops2
from cLoops2.ds import Peak, Loop, XY
from cLoops2.io import ixy2pet, parseIxy, writeNewJson, parseBed2Peaks, parseTxt2Loops


def stichRegions(rs, gap=1):
    """
    Stich 1D regions with specified gap size. 
    """
    cov = set()
    for r in rs:
        cov.update(range(r[0], r[1] + 1))
    cov = list(cov)
    cov.sort()
    nrs = []
    i = 0
    while i < len(cov) - 1:
        for j in range(i + 1, len(cov)):
            if cov[j] - cov[j - 1] > gap:
                break
            else:
                continue
        start = cov[i]
        end = cov[j - 1]
        nrs.append([start, end])
        i = j  #update search start
    return nrs


def filterPETs(rs, key, predir, fixy, iv=False):
    """
    Filter PETs, only keep those located at regions.
    """
    print("%s\t Filtering PETs of %s with %s regions." %
          (datetime.now(), key, len(rs)))
    key2, mat = parseIxy(fixy)
    xy = XY(mat[:,0],mat[:,1])
    rids = set()
    for r in tqdm(rs):
        r = xy.queryPeak(r[0], r[1])
        rids.update(r)
    rids = list(rids)
    if len(rids) == 0:
        return
    if iv:
        aids = set(np.arange(mat.shape[0]))
        rids = aids.difference(rids)
    mat = mat[rids, ]
    foixy = predir + "/" + "-".join(key2) + ".ixy"
    joblib.dump(mat, foixy)


def filterPETsByPeaks(predir, fbed, outdir, cpu=1, iv=False, gap=1):
    """
    Filter PETs according to peaks form .bed file. 
    """
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    peaks = parseBed2Peaks(fbed)
    npeaks = {}
    for key in peaks:
        if key not in meta["data"]["cis"]:
            continue
        rs = [[p.start, p.end] for p in peaks[key]]
        npeaks[key] = stichRegions(rs, gap)
    Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(filterPETs)(
        npeaks[key],
        key,
        outdir,
        meta["data"]["cis"][key]["ixy"],
        iv,
    ) for key in npeaks.keys())
    writeNewJson(outdir)


def _filterPETsByLoops(loops,key,predir,fixy,iv=False):
    """
    Filter PETs by loops
    """
    print("%s\t Filtering PETs of %s with %s regions." %
          (datetime.now(), key, len(loops)))
    key2, mat = parseIxy(fixy)
    xy = XY(mat[:,0],mat[:,1]) 
    rids = set()
    for loop in tqdm(loops):
        a, b,r = xy.queryLoop(loop.x_start, loop.x_end, loop.y_start, loop.y_end)
        rids.update(r)
    rids = list(rids)
    if len(rids) == 0:
        return
    if iv:
        aids = set(np.arange(mat.shape[0]))
        rids = aids.difference(rids)
    mat = mat[rids, ]
    foixy = predir + "/" + "-".join(key2) + ".ixy"
    joblib.dump(mat, foixy)



def filterPETsByLoops(predir, floop, outdir, cpu=1, iv=False, gap=1,both=False):
    """
    Filter PETs according to loops from _loop.txt file. 
    """
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    loops = parseTxt2Loops(floop)
    if both:
        #filter PETs, only keep those both ends overlapped with target loop anchors
        Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(_filterPETsByLoops)(
            loops[key],
            key,
            outdir,
            meta["data"]["cis"][key]["ixy"],
            iv,
        ) for key in loops.keys())
    else:
        #filter PETs, keep those any end overlapped with target loop anchors
        npeaks = {}
        for key in loops:
            if key not in meta["data"]["cis"]:
                continue
            rs = []
            for loop in loops[key]:
                rs.append([loop.x_start, loop.x_end])
                rs.append([loop.y_start, loop.y_end])
            npeaks[key] = stichRegions(rs, gap)
        Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(filterPETs)(
            npeaks[key],
            key,
            outdir,
            meta["data"]["cis"][key]["ixy"],
            iv,
        ) for key in npeaks.keys())
    writeNewJson(outdir)


def _filterPETsBySingletons(f, outdir, binSize):
    """
    @param f:str .ixy file
    @param outir: str,
    @param bs: int, binSize
    """
    key, mat = parseIxy(f)
    print("%s\t Filtering %s singleton PETs in contact matrix bins." %
          (datetime.now(), key))
    minC = np.min(mat)
    ss = {}
    i = 0
    for x,y in tqdm(mat):
        x = int((x - minC) / binSize)
        y = int((y - minC) / binSize)
        if x not in ss:
            ss[x] = {}
        if y not in ss[x]:
            ss[x][y] = []
        ss[x][y].append(i)
        i += 1
    rs = []
    for nx in ss.keys():
        for ny in ss[nx].keys():
            if len(ss[nx][ny]) > 1:
                rs.extend(ss[nx][ny])
    if len(rs) > 0:
        mat = mat[rs, ]
        foixy = outdir + "/" + "-".join(key) + ".ixy"
        joblib.dump(mat, foixy)


def filterPETsBySingletons(predir, outdir, bs, cpu=1):
    """
    Filter singleton PETs in contact matrix bins.
    """
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(_filterPETsBySingletons)(
        meta["data"]["cis"][key]["ixy"],
        outdir,
        bs,
    ) for key in meta["data"]["cis"].keys())
    writeNewJson(outdir)


def _getNearbyGrids(Gs, cell):
    x, y = cell[0], cell[1]
    keys = [(x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y), (x - 1, y - 1),
            (x - 1, y + 1), (x + 1, y - 1), (x + 1, y + 1)]
    ncells = []
    for key in keys:
        if key in Gs:
            ncells.append(key)
    return ncells


def _filterPETsByKNNs(f, outdir,eps,minPts):
    """
    @param f:str .ixy file
    @param outir: str,
    """
    key, mat = parseIxy(f)
    print("%s\t Filtering %s PETs based on KNN." % (datetime.now(), key))
    
    #build grids
    minX, minY = np.min(mat[:,0]),np.min(mat[:,1])
    Gs = {}
    ps = {}
    for i,(x,y) in enumerate(mat):
        nx = int((x - minX) / eps) + 1
        ny = int((y - minY) / eps) + 1
        Gs.setdefault((nx, ny), []).append(i)
        #last elements marks the class, initially -1 as noise
        ps[ i ] = [x, y, nx, ny, -1]

    #Grid index with all neighbor points.
    Gs2 = {}
    for cell in Gs.keys():
        nps = []
        nps.extend(Gs[cell])
        for cellj in _getNearbyGrids(Gs,cell):
            nps.extend(Gs[cellj])
        Gs2[cell] = nps
    
    #remove noise 
    #: noise cells without neighbors
    tode = set()
    #: noise cells with neighbors
    tode2 = set()
    for cell in Gs.keys():
        if len(Gs2[cell]) < minPts:
            tode2.add(cell)
    #KNN to noise cells with neighbors
    for cell in tode2:
        cells = _getNearbyGrids(Gs,cell)
        ncells = set(cells) & tode2
        #all neighbor cells are noise
        if len(cells) == len(ncells):
            tode.add(cell)
    for cell in tode:
        for p in Gs[cell]:
            del ps[p]
        del Gs[cell]
    
    nps = []
    for cell in Gs.keys():
        nps.extend( Gs[cell] )

    if len(nps) > 0:
        mat = mat[nps,]
        foixy = outdir + "/" + "-".join(key) + ".ixy"
        joblib.dump(mat, foixy)


def filterPETsByKNNs(predir, outdir, eps=1000, minPts=5, cpu=1):
    """
    Filter PETs based on blockDBSCAN noise-removing processing.
    """
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(_filterPETsByKNNs)(
        meta["data"]["cis"][key]["ixy"],
        outdir,
        eps,
        minPts,
    ) for key in meta["data"]["cis"].keys())
    writeNewJson(outdir)


def _samplePETs(f, outdir, r):
    """
    @param f: str .ixy file
    @param outir: str, output directory
    @param r: float, sampling ratio
    """
    key, mat = parseIxy(f)
    print("%s\t Sampling PETs for %s." % (datetime.now(), key))
    tr = int(round(mat.shape[0] * r))
    if tr < 1:
        return
    if tr < mat.shape[0]:
        rs = np.random.choice(mat.shape[0], tr,replace=False)
    else:
        rs = np.random.choice(mat.shape[0], tr,replace=True)
    mat = mat[rs, ]
    #update PET id
    foixy = outdir + "/" + "-".join(key) + ".ixy"
    joblib.dump(mat, foixy)


def samplePETs(predir, outdir, tot, cpu=1):
    """
    Sample PETs to target library size.
    """
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    uni = meta["Unique PETs"]
    r = tot / float(uni)
    Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(_samplePETs)(
        meta["data"]["cis"][key]["ixy"],
        outdir,
        r,
    ) for key in meta["data"]["cis"].keys())
    writeNewJson(outdir)
    """
    Parallel(n_jobs=cpu,backend="multiprocessing")(delayed(_samplePETs)(
        meta["data"]["trans"][key]["ixy"],
        outdir,
        r,
    ) for key in meta["data"]["trans"].keys() )
    """
