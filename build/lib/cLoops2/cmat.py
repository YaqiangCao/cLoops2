#!/usr/bin/env python
#--coding:utf-8--
"""
cmat.py
cLoops2 contact matrix and 1D signal pileup related functions.
2019-12-17: 1D track support added
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#3rd
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.sparse import csr_matrix

#cLoops2
from cLoops2.ds import XY
from cLoops2.io import parseIxy


def xy2dict(mat, s=-1, e=-1, r=5000):
    """
    Convert the coordinates to contact matrix with specified resolution.
    @param mat: np.array, [[x,y]]
    @param s: start site to build the contact matrix, if -1, infer from the data min , assigned can help to maintain same shape
    @param e: end site to build the contact matrix, if -1, infer from the data max , assigned can help to maintain same shape
    @param r: resolution for the contact matrix
    """
    nmat = {}
    if s == -1:
        s = np.min(mat)
    if e == -1:
        e = np.max(mat)
    for x, y in mat:
        nx = int((x - s) / r)
        ny = int((y - s) / r)
        if nx not in nmat:
            nmat[nx] = {}
        if ny not in nmat[nx]:
            nmat[nx][ny] = 0
        nmat[nx][ny] += 1
    #add the max value to maintain the same shape of matrix across data
    nx = int((e - s) / r)
    ny = int((e - s) / r)
    if nx not in nmat or ny not in nmat[nx]:
        nmat[nx] = {ny: 0}
    return nmat


def ixy2dict(f, r=5000, cut=0):
    """
    Convert the .ixy file to contact matrix with specified resolution.
    @param f: .ixy file
    @param r: resolution for the contact matrix
    @param cut: distance cutoff to filter PETs
    """
    chrom, mat = parseIxy(f, cut=cut)
    return xy2dict(mat, r=r)


def dict2mat(nmat):
    """
    Conver contact matrix in dict to sparse matrix.
    @return scipy.sparse import csr_matrix
    """
    data, row, col = [], [], []
    for nx in nmat.keys():
        for ny in nmat[nx].keys():
            #create the symetric matrix
            data.append(nmat[nx][ny])
            row.append(nx)
            col.append(ny)
            data.append(nmat[nx][ny])
            row.append(ny)
            col.append(nx)
    cmat = csr_matrix((data, (row, col)))
    return cmat


def getObsMat(xy, start, end, r):
    """
    Get the observed interaction contact matrix.
    xy is [[x,y]]
    r is resolution
    """
    ps = np.where(xy[:, 0] >= start)[0]
    xy = xy[ps, ]
    ps = np.where(xy[:, 1] <= end)[0]
    xy = xy[ps, ]
    mat = xy2dict(xy, s=start, e=end, r=r)
    mat = dict2mat(mat)
    mat = mat.toarray()
    return mat


def getExpMat(xy, shape, start, end, r, repeats=5):
    """
    Get the expected interaction contact matrix.
    xy is [[x,y]]
    shape is () shape from the observed matrix.
    r is resolution
    """
    mat = []
    i = 0
    while i < repeats:
        a = xy[:, 0]
        b = xy[:, 1]
        np.random.shuffle(a)
        np.random.shuffle(b)
        xy[:, 0] = a
        xy[:, 1] = b
        s = b-a
        s = np.where( s > 0)[0]
        nxy = xy[s,] 
        nmat = getObsMat(nxy, start, end, r)
        if nmat.shape == shape:
            mat.append(nmat)
            i += 1
    mat = np.array(mat)
    return mat.mean(axis=0)


def get1DSig(xy, start, end, ext=50):
    """
    Get the overlayed 1D signal
    @param xy, cLoops2.ds.XY object
    @param start: int, start coordinate
    @param end: int, end coordinate
    @param ext: int, extention of each tag
    """
    ss = np.zeros(end - start)
    l_idx = np.searchsorted(xy.xs, start, side="left")
    r_idx = np.searchsorted(xy.xs, end, side="right")
    for i in range(l_idx, r_idx):
        x = xy.xs[i]
        pa = max(0, x - start - ext)
        pb = min(max(0, x - start + ext), end-start) #fix max
        ss[pa:pb] += 1
    l_idx = np.searchsorted(xy.ys, start, side="left")
    r_idx = np.searchsorted(xy.ys, end, side="right")
    for i in range(l_idx, r_idx):
        y = xy.ys[i]
        pa = max(0, y - start - ext)
        pb = min(max(0, y - start + ext), end)
        ss[pa:pb] += 1
    return ss


def getBinMean(s, bins=100):
    """
    Get the mean of bins for a array.
    @param s: np.array
    @param bins: int, how many bins as converted
    """
    width = int(len(s) / bins)
    ns = s[:bins * width].reshape(-1, width).mean(axis=1)
    return ns


def get1DSigMat(xy, rs, ext=20, bins=100, skipZeros=False):
    """
    Get the 1D signal matrix for a set of regions.
    @param xy is XY object
    ext is extend of reads from the center
    bins is the final array size for a record
    return a pd.Dataframe, row is regions/peaks, columns is j
    """
    ds = {}
    print("Get 1D signal for %s regions" % (len(rs)))
    for r in tqdm(rs):
        #for r in rs:
        if r[2] - r[1] < bins:
            continue
        s = get1DSig(xy, int(r[1]), int(r[2]), ext=ext)
        if skipZeros and np.sum(s) == 0:
            continue
        ns = getBinMean(s, bins=bins)
        if len(r) >= 4:
            rid = "|".join(list(map(str, r[:4])))
        else:
            rid = "|".join(list(map(str, r[:3])))
        ds[rid] = ns
    if len(ds) == 0:
        return None
    else:
        ds = pd.DataFrame(ds).T
        return ds



def getVirtual4CSig(xy,start,end,viewStart,viewEnd,ext=20):
    """
    Get the virtual 4C signal for a region and a view point.
    @param xy, cLoops2.ds.XY object
    @param start: int, start coordinate
    @param end: int, end coordinate
    @param viewStart: int, start coordinate for view point
    @param viewEnd: int, end coordinate for view point
    """
    aps = xy.queryPeak(start,end)
    bps = xy.queryPeak(viewStart,viewEnd)
    cps = xy.queryPeakBoth(viewStart,viewEnd)
    ps = aps.intersection( bps ).difference( cps )
    ss = np.zeros(end - start)
    for i in ps:
        x = xy.mat[i,0]
        y = xy.mat[i,1]
        pa = max(0, x - start - ext)
        pb = min(max(0, x - start + ext), end-start) #fix max
        ss[pa:pb] += 1
        pa = max(0, y - start - ext)
        pb = min(max(0, y - start + ext), end)
        ss[pa:pb] += 1
    return ss



