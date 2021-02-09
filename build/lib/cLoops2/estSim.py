#!/usr/bin/env python
#--coding:utf-8 --
"""
cLoops2:estSim
cLoops2 similarity estimation of 3D genome interactions.
The HiC-Spector method is too time cosuming for measuring the similarity, so we turn to other methods, especially for Trac-looping data. 
To get the laplacian normalized matrix, can be easily through mat = scipy.aparse.csgraph.laplacian(mat,normed=False), and the eigenvectors can be obtained through scipy.sparse.linalg.eigs
2019-09-24: PLS fails as decomposition with all nan, CCA fails as too time costly. PCA works well for trac-looping.
2019-09-24: finished. 
2020-06-25: extend to multiple data sets
2020-08-19: add PETs cutoff
"""

__date__ = "2019-09-17"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import warnings
warnings.filterwarnings("ignore")
import os
import json
import random
from glob import glob

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn.decomposition import PCA

#cLoops2
from cLoops2.io import parseIxy
from cLoops2.cmat import xy2dict, dict2mat, getObsMat, getExpMat
from cLoops2.plot import plotCorrScatterPCC, plotCorrScatterPCA, plotCorrHeatmap
from cLoops2.settings import *


def pre(dirs, cut, mcut):
    """
    Prepare samples.
    """
    #prepare files
    chroms = {}
    #for all samples find the
    for dir in dirs:
        sample = dir.split("/")[-1]
        metaf = dir + "/petMeta.json"
        meta = json.loads(open(metaf).read())
        for k, v in meta["data"]["cis"].items():
            if k not in chroms:
                chroms[k] = {"samples": {}, "start": np.inf, "end": -1}
            fixy = v["ixy"]
            chroms[k]["samples"][sample] = fixy
            key, mat = parseIxy(fixy, cut=cut, mcut=mcut)
            s = np.min(mat)
            e = np.max(mat)
            if s < chroms[k]["start"]:
                chroms[k]["start"] = s
            if e > chroms[k]["end"]:
                chroms[k]["end"] = e
    #just in case some samples do not have the chrom
    for dir in dirs:
        sample = dir.split("/")[-1]
        metaf = dir + "/petMeta.json"
        meta = json.loads(open(metaf).read())
        for k in chroms.keys():
            if k not in meta["data"]["cis"]:
                chroms[k]["samples"][sample] = None
    #chrom, key is chr1-chr1, values is {sampleid:fixy,start:start,end:end}
    return chroms


def getFlattenVector(key, samples, start, end, r=5000, cut=0, mcut=-1):
    """
    Flatten the contact matrix as vectors for multiple samples the same chrom-chrom.
    """
    print("Getting the data of %s, from %s samples." % (key, len(samples)))
    #all possible locations
    xykeys = set()
    #contact matrix and flatten vector
    obsMat, obsVec = {}, {}
    #get all contact matrix values
    for sample, fixy in tqdm(samples.items()):
        if fixy is not None:
            nkey, mat = parseIxy(fixy, cut=cut, mcut=mcut)
        else:
            mat = np.array([])
        mat = xy2dict(mat, s=start, e=end, r=r)
        for nx in mat.keys():
            for ny in mat[nx].keys():
                xykeys.add((nx, ny))
        obsMat[sample] = mat
        obsVec[sample] = []
    #flatten them
    print("Combing the data of %s." % key)
    for key in xykeys:
        #observation
        for sample in obsMat.keys():
            if key[0] in obsMat[sample] and key[1] in obsMat[sample][key[0]]:
                obsVec[sample].append(obsMat[sample][key[0]][key[1]])
            else:
                obsVec[sample].append(0)
    return obsVec


def comparePCC(dirs, fout, r=5000, cut=0, mcut=-1, cpu=1, pcut=2, plot=False):
    """
    Caculating the PCC of flatten vectors from contact matrix for similarity.
    """
    chroms = pre(dirs, cut, mcut)
    print("Estimating PCC for %s." % (",".join(dirs)))
    ds = Parallel(n_jobs=cpu,
                  backend="multiprocessing")(delayed(getFlattenVector)(
                      k,
                      v["samples"],
                      v["start"],
                      v["end"],
                      r=r,
                      cut=cut,
                      mcut=mcut,
                  ) for k, v in chroms.items())
    obsMat = {}
    for d in ds:
        for s, v in d.items():
            if s not in obsMat:
                obsMat[s] = []
            obsMat[s].extend(v)
    del ds
    obsMat = pd.DataFrame(obsMat)
    ns = []
    for t in obsMat.itertuples():
        s = np.array(t[1:])
        s = s[s > pcut]
        if len(s) == 0:
            ns.append(t[0])
    obsMat = obsMat.drop(ns)
    obsMat.to_csv(fout + "_PCC_obsVectors.txt", sep="\t", index_label="binId")
    obsCorr = obsMat.corr()
    obsCorr.to_csv(fout + "_PCC.txt", sep="\t", index_label="sample")
    if plot:
        if len(dirs) == 2:
            #density scatter plot
            plotCorrScatterPCC(obsMat, fout + "_PCC.pdf")
        else:
            #heatmap
            plotCorrHeatmap(obsCorr, fout + "_PCC.pdf")


def getPCAFlattenVector(key,
                        samples,
                        start,
                        end,
                        r=5000,
                        n_comp=2,
                        cut=0,
                        mcut=-1):
    """
    Flatten the contact matrix as vectors for multiple samples the same chrom-chrom.
    @param s: int, start site for construction of contact matrix
    @param e: int, end site for construction of contact matrix
    @param r: int, resolution 
    @param n_comp: top n components
    """
    print("Getting the data of %s, from %s samples." % (key, len(samples)))
    amat = None
    #get all contact matrix values, do not keep in case of out-of-memory
    for i, (sample, fixy) in enumerate(tqdm(samples.items())):
        if fixy is not None:
            nkey, mat = parseIxy(fixy, cut=cut, mcut=mcut)
        else:
            mat = np.array([])
        mat = xy2dict(mat, s=start, e=end, r=r)
        mat = dict2mat(mat)
        mat = mat.todense()
        if i == 0:
            amat = mat
        else:
            amat = amat + mat
    print("Getting the projected PCs of %s." % key)
    #average contact matrix and PCA projection
    amat = amat / len(samples)
    pca = PCA(n_components=n_comp)
    pca.fit(amat)
    del amat
    #PCA flatten vector
    obsVec = {}
    for sample, fixy in tqdm(samples.items()):
        key, mat = parseIxy(fixy, cut=cut, mcut=mcut)
        mat = xy2dict(mat, s=start, e=end, r=r)
        mat = dict2mat(mat)
        mat = mat.todense()
        mat = pca.transform(mat)
        obsVec[sample] = mat.flatten()
    return obsVec


def comparePCA(dirs, fout, r=5000, cut=0, mcut=-1, cpu=1, n_comp=2,
               plot=False):
    """
    Caculating the PCC of flatten vectors from contact matrix PCA for similarity.
    """
    chroms = pre(dirs, cut, mcut)
    print("Estimating PCC for %s." % (",".join(dirs)))
    ds = Parallel(n_jobs=cpu,
                  backend="multiprocessing")(delayed(getPCAFlattenVector)(
                      k,
                      v["samples"],
                      v["start"],
                      v["end"],
                      r=r,
                      n_comp=n_comp,
                      cut=cut,
                      mcut=mcut,
                  ) for k, v in chroms.items())
    obsMat = {}
    for d in ds:
        for s, v in d.items():
            if s not in obsMat:
                obsMat[s] = []
            obsMat[s].extend(v)
    del ds
    obsMat = pd.DataFrame(obsMat)
    obsMat.to_csv(fout + "_PCA_obsVectors.txt", sep="\t", index_label="binId")
    obsCorr = obsMat.corr()
    obsCorr.to_csv(fout + "_PCA.txt", sep="\t", index_label="sample")
    if plot:
        if len(dirs) == 2:
            #density scatter plot
            plotCorrScatterPCA(obsMat, fout + "_PCA.pdf")
        else:
            #heatmap
            plotCorrHeatmap(obsCorr, fout + "_PCA.pdf")


def estSim(dirs,
           fout,
           bs=5000,
           method="pcc",
           cut=0,
           mcut=-1,
           cpu=1,
           n_comp=2,
           pcut=2,
           plot=False):
    """
    Estimate interaction similarities.
    @param dirs: list of str, cLoops2 pre generated data directories
    @param fout: str, output prefix
    @param bs: int, binSize for contact matrix 
    @param method: str, options are pcc and pca
    @param cut: int, >cut distance PETs will be used 
    @param mcut: int, <mcut distance PETs will be used 
    @param cpu: int, cpu numbers to run jobs
    @param n_comp: int, first n_comp components were used for PCA embeding 
    @param pcut: int,for a bin, if all samples in that bin pets <= pcut, remove the bin
    @param plot: bool, whether to plot the result
    """
    if method == "pcc":
        comparePCC(
            dirs,
            fout,
            r=bs,
            cut=cut,
            mcut=mcut,
            cpu=cpu,
            pcut=pcut,
            plot=plot,
        )
    elif method == "pca":
        comparePCA(
            dirs,
            fout,
            r=bs,
            cut=cut,
            mcut=mcut,
            cpu=cpu,
            n_comp=n_comp,
            plot=plot,
        )
    else:
        print("ERROR! The selected method not implemented! Return.")
