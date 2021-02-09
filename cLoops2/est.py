#!/usr/bin/env python
#--coding:utf-8--
"""
est.py
Estimate key parameters/cutoffs/models for cLoops2 specific models.
2020-12-18: TheilSenRegressor if fit_intercept is set to False, then predict will have problem.
"""

__author__ = "CAO Yaqiang"
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os
import json

#3rd
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn import linear_model
from joblib import Parallel, delayed
from sklearn.mixture import GaussianMixture as GMM
from sklearn.neighbors import NearestNeighbors as NN

#sklearn
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import mean_absolute_error as MAE

#keras
"""
os.environ["TF_CPP_MIN_LOG_LEVEL"] = '3'
import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras import metrics
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import load_model, Sequential, Model
from tensorflow.keras.callbacks import ModelCheckpoint, ReduceLROnPlateau, EarlyStopping
from tensorflow.keras.layers import Flatten, Dense, Dropout, BatchNormalization, Activation
"""

#warning settings
import warnings
warnings.filterwarnings("ignore")

#cLoops2
from cLoops2.io import parseIxy
from cLoops2.plot import plotEstRes, plotEstSat
from cLoops2.settings import *


#### Estimation of eps based on expected domain knowedge that 3D interaction data may have local peaks, limited to use intra-chromosomal PETs
def getXyDis(f, cut=0, mcut=-1):
    """
    Get the distance between PETs, X-Y for a .ixy file.
    """
    key, mat = parseIxy(f, cut, mcut)
    j = mat.shape[0]
    if j < 2:
        return None
    dis = np.abs(mat[:, 1] - mat[:, 0])
    dis = np.array(dis)
    dis = dis[~np.isnan(dis)]
    dis = dis[~np.isinf(dis)]
    dis = dis[dis > 0]
    return dis


def getGmmLabelsEps(dis, n_components=2):
    """
    Estimate Gaussian Mixture Models.
    @param dis: np.array of distance between (X,Y), log2 transformed
    @param n_components: int, expected components number for the mixture model, as expect 2, self-ligation and inter-ligation two classes
    @return: ps, numpy.array, indicate the class of point; eps, estimated eps 
    """
    #dis = np.log2(dis)
    ds = dis.reshape(dis.shape[0], 1)
    cvs = ['spherical', 'tied', 'diag', 'full']  #covariance types
    lowest_bic = np.infty
    clf = None
    for cv in cvs:
        gmm = GMM(n_components=n_components,
                  covariance_type=cv,
                  random_state=123,
                  weights_init=[0.5, 0.5])
        gmm.fit(ds)
        bic = gmm.bic(ds)
        if bic < lowest_bic:
            lowest_bic = bic
            clf = gmm
    #print(clf.means_)
    eps = 2 * int(2**min(clf.means_))
    ps = np.array(clf.predict(ds))  #predicted labels of PETs from clf
    return ps, eps


### Estimation of eps based on k-distance plot, no limitations, maybe not work well
def getKDis(f, k=5, cut=0, cpu=1):
    """
    Get the k-nearest neighbor k-distance based on sklearn.neighbors.NearestNeighbors.
    @param f: str,.ixy file
    @param k: int, the k neighbor
    @param cpu: int, the cpu numbers to run the job
    @param cut: int, the distance cutoff to filter PETs.
    """
    key, mat = parseIxy(f, cut)
    j = mat.shape[0]
    if j < 2:
        return None
    nn = NN(n_neighbors=k,
            algorithm='ball_tree',
            metric="cityblock",
            n_jobs=cpu).fit(mat)
    dis, indices = nn.kneighbors(mat)
    dis = dis[:, -1]
    return dis


def getKDisKneeEps(dis, S=10):
    """
    Estimate the knee point from the ascending k-distance plot
    @param dis: numpy.array, ascending sorted log2 k-dis 
    @param S: int, the dis data bin to how many slices
    @return: the position of the knee and the estimated eps
    """
    step = int(len(dis) / S)
    maxi, maxfc = -1, 0
    for i in range(2 * step, len(dis) - step):  #skip the fisrt knee
        fc = np.mean(dis[i + 1:i + step]) / np.mean(
            dis[i - step:i])  #for a point, check the upstream/downstream
        if fc > maxfc:
            maxi = i
            maxfc = fc
    eps = dis[maxi]
    return maxi, eps


### Estimation of distance cutoff for intra-chromosomal inter-ligaiton vs. self-ligation PETs
def solve(m1, std1, m2, std2):
    """
    Finding the intersection point of two gaussian curves
    From: https://stackoverflow.com/questions/22579434/python-finding-the-intersection-point-of-two-gaussian-curves
    """
    a = 1 / (2 * std1**2) - 1 / (2 * std2**2)
    b = m2 / (std2**2) - m1 / (std1**2)
    c = m1**2 / (2 * std1**2) - m2**2 / (2 * std2**2) - np.log(std2 / std1)
    r = np.roots([a, b, c])
    r = r[r > 0]
    #r = r[r > 10]  #just used to avoid too small cutoffs, may have bias for Trac-looping
    return min(r)


def estIntraCut(di, ds, log=True):
    """
    Estimation of distance cutoff for inter-ligation and self-ligation pets.
    @param di: list,distance for inter-ligation cluster pets
    @param ds: list,distance for self-ligation cluster pets
    """
    di = np.abs(np.array(di))
    ds = np.abs(np.array(ds))
    di = di[~np.isnan(di)]
    ds = ds[~np.isnan(ds)]
    di = di[di > 0]
    ds = ds[ds > 0]
    if log:
        di = np.log2(di)
        ds = np.log2(ds)
    cut = [np.median(ds) + 3 * ds.std()]
    cut.append(solve(di.mean(), di.std(), ds.mean(), ds.std()))
    cut = [ c for c in cut if c < di.mean()]
    cut = max(cut)
    #cut = min(cut)
    if log:
        cut = int(2**cut)
    return cut


## Estimation of distance and interaction density
def getObsDisFreq(mat, binSize=1000):
    """
    Get the relation between genomic distance with interactions using bin size based on Numpy
    """
    minC = np.min(mat)
    a = (mat[:, 0] - minC) / binSize
    b = (mat[:, 1] - minC) / binSize
    a = a.astype(int)
    b = b.astype(int)
    c = np.abs(a - b) * binSize
    c = c.astype(int)
    sso = Counter(c)
    return sso


def preObs(fixy, cut=0, binSize=1000):
    """
    Process observed data.
    """
    chrom, mat = parseIxy(fixy, cut=cut)
    sso = getObsDisFreq(mat, binSize=binSize)
    return sso


def updateFreq(ssa, ssb):
    """
    Update the frequency dict
    """
    for k, v in ssb.items():
        if k not in ssa:
            ssa[k] = v
        else:
            ssa[k] = ssa[k] + v
    return ssa


def estLr(meta, cpu=1, plot=False, fout=None):
    """
    Estimation the linear relationship between genomic distance and interaction frequency.
    """
    ds = Parallel(n_jobs=cpu, backend="multiprocessing")(
        delayed(preObs)(meta["data"]["cis"][key]["ixy"], cut=0, binSize=10000)
        for key in meta["data"]["cis"].keys())
    rs = {}
    for d in ds:
        rs = updateFreq(rs, d)
    del rs[0]  #remove distance 0
    x = np.log10(np.array(list(rs.keys())))
    y = np.array(list(rs.values()))
    y = y / y.sum()
    y = np.log10(y)
    rs = pd.Series(y, index=x)
    #only using reads within 10M get better fitting line
    nx = x[x > 0]
    nx = nx[nx < 7]
    ny = rs[nx].values
    nnx = np.array([[t] for t in nx])
    lr = linear_model.LinearRegression(n_jobs=cpu, fit_intercept=True)
    lr.fit(nnx, ny)
    if plot and fout is not None:
        fig, ax = pylab.subplots()
        ax.scatter(x, y, color=colors[0], s=2, label="observed")
        ax.scatter(nx,
                   lr.predict(nnx),
                   color=colors[1],
                   s=2,
                   label="linear fit")
        ax.legend()
        ax.set_xlabel("Genomic distance, log10(bp)")
        ax.set_ylabel("Normalized interaction frequency,log10")
        ax.set_title("y=%sx+%s" % (lr.coef_[0], lr.intercept_))
        pylab.tight_layout()
        pylab.savefig(fout + "_disFreqFit.pdf")
    #the lr could be used as
    #n = int(N * 10**(lr.predict([[np.log10(loop.distance)]])[0]))  #the potential PETs that share the similar distance
    return lr


## Estimation of contact matrix resolution
def getBinPETs(f, binSize, cut=0, mcut=-1):
    """
    Get the number of PETs in bins.
    @param f: str, .ixy file
    @param binSize:int, contact matrix bin size
    """
    chrom, mat = parseIxy(f, cut, mcut)
    print("Get signals from", chrom)
    minC = np.min(mat)
    ss = {}
    for x, y in tqdm(mat):
        x = int((x - minC) / binSize)
        y = int((y - minC) / binSize)
        if x == y:  #important to remove self-ligation
            continue
        if x not in ss:
            ss[x] = {}
        if y not in ss[x]:
            ss[x][y] = 0
        ss[x][y] += 1
    sso = []
    for x in ss.keys():
        for y in ss[x].keys():
            sso.append(ss[x][y])
    return np.array(sso)


def getCumBins(ds, bins=100):
    """
    Furthur bin the signal in contact matrix into bins, only care of the cumutative trend.
    """
    nn = []
    step = int(len(ds) / bins)
    for i in range(0, len(ds), step):
        if i + step > len(ds):
            break
        nn.append(ds[i:i + step].sum())
    nn = np.array(nn)
    nn = np.cumsum(nn) / float(nn.sum()) * 100
    nn = pd.Series(nn, index=list(range(bins)))
    return nn


def estRes(predir,
           fnOut,
           logger,
           bs=[25000, 5000, 1000],
           cpu=1,
           cut=0,
           mcut=-1):
    """
    Estimation of reasonable contact matrix resolution based on 2D signal enrichment.
    """
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    fs = []
    for key in meta["data"]["cis"]:
        fs.append(meta["data"]["cis"][key]["ixy"])
    cumBins = []
    singletonRatios = []
    PETsRatio = []
    data = {}
    for binSize in bs:
        print("Get the signal distribution with resolution of %s" % binSize)
        #ds = Parallel(n_jobs=cpu, backend='threading')(delayed(getBinPETs)(
        ds = Parallel(n_jobs=cpu, backend='multiprocessing')(
            delayed(getBinPETs)(f, binSize, cut=cut, mcut=mcut) for f in fs)
        ds = np.concatenate(ds)
        #sort the data and get how many are singletons
        ds = np.sort(ds)
        p = np.where(ds[ds <= 1])[0]
        t = len(p) / np.sum(ds) * 100
        p = np.max(p)
        r = p / float(len(ds)) * 100
        ss = getCumBins(ds)
        cumBins.append(ss)
        singletonRatios.append(r)
        PETsRatio.append(t)
        logger.info(
            "resolution %s:  %.2f%% contact matrix bins only contain singleton PET, singleton PETs ratio %.2f%%."
            % (binSize, r, t))
        data[ binSize ] = ss
    data = pd.DataFrame(data)
    data.to_csv(fnOut+"_estRes.txt",sep="\t",index_label="percentage of bins")
    plotEstRes(bs, cumBins, singletonRatios, PETsRatio, fnOut)


## Estimation of scaling factor between two samples
def estSf(cs, ts, cpu=1):
    """
    Estimate the scaling factor using linear regression or the targetTotal/controlTotal
    @param controlCounts: list of int, paired with targetCounts
    @param targetCounts: list of int, paired with controlCounts
    @param cpu: cpu number to run the fitting jobs
    """
    x = [[c] for c in cs]
    x = np.array(x)
    ts = np.array(ts)
    x_train, x_vali, y_train, y_vali = train_test_split(x, ts, test_size=0.1)
    lra = linear_model.LinearRegression(n_jobs=cpu, fit_intercept=False)
    lrb = linear_model.TheilSenRegressor(random_state=123,
                                         n_jobs=cpu,
                                         fit_intercept=False)
    lra.fit(x_train, y_train)
    lrb.fit(x_train, y_train)
    sfa = lra.coef_[0]
    sfb = lrb.coef_
    xv = np.array([ t[0] for t in x_vali])
    ypa = xv * sfa 
    ypb = xv * sfb
    #ypa = lra.predict(x_vali)
    #ypb = lrb.predict(x_vali)
    maea = MAE(y_vali, ypa)
    maeb = MAE(y_vali, ypb)
    if maea < maeb:
        sf = lra.coef_[0]
        lr = lra
    else:
        sf = lrb.coef_
        lr = lrb
    """
    #plot function used for debug
    fig, ax = pylab.subplots()
    xs = [ t[0] for t in x_vali]
    ax.scatter( xs, y_vali,label="validataion",s=1,color=colors[0] )
    ax.scatter( xs, lr.predict(x_vali),label="fitted,sf:%.3f"%(sf),s=1,color=colors[1])
    ax.legend()
    pylab.savefig(fout+"_sf.pdf")
    """
    return sf


def estSfMANorm(cs, ts):
    """
    Estimate the linear fitting of background for target and control based on 
    the assumption of MANorm2. 
    @param controlCounts: list of int, paired with targetCounts, 
        interaction density per kb, already log2 transformed
    @param targetCounts: list of int, paired with controlCounts
    """
    beta = np.std(cs) / np.std(ts)
    alpha = np.mean(cs) - beta * np.mean(ts)
    return [beta, alpha]


def dEstSf(cs, ts, fout, cpu=1):
    """
    Estimate scaling factor based on deep-learning. 
    """
    #x = [[c] for c in cs]
    x_train, x_vali, y_train, y_vali = train_test_split(cs, ts, test_size=0.1)
    model = Sequential()
    model.add(Dense(32, input_dim=1))
    model.add(BatchNormalization())
    model.add(Activation(activation="sigmoid"))
    model.add(Dense(16, activation="relu"))
    model.add(Dropout(0.1))
    model.add(Dense(1, activation="relu"))
    model.compile(loss='mae', optimizer=Adam(1e-3), metrics=['mse', 'mae'])
    reduce_lr = ReduceLROnPlateau(monitor="val_loss", patience=10)
    early_stop = EarlyStopping(monitor='val_loss', patience=20)
    callbacks = [reduce_lr, early_stop]

    hist = model.fit(x=x_train,
                     y=y_train,
                     callbacks=callbacks,
                     epochs=100,
                     shuffle=True,
                     validation_data=(x_vali, y_vali))

    yps = model.predict(x_vali)
    yps = [t[0] for t in yps]
    """
    fig, ax = pylab.subplots()
    ax.scatter(x_vali,y_vali,label="validation",s=1,color=colors[0] )
    ax.scatter(x_vali,yps,label="fitted_validation",s=1,color=colors[1] )
    ax.legend()
    pylab.savefig(fout+"_dsf.pdf")
    """
    return model


## Estimation of sequencing saturation
def getSampBins(f, binSize, cut=0, mcut=-1, repeats=3, tol=5):
    """
    Get the number of detect bins from sub-sampling.
    @param f: str, .ixy file
    @param binSize:int, contact matrix bin size
    """
    chrom, mat = parseIxy(f, cut, mcut)
    if mat.shape[0] == 0:
        return None
    print("Get signals from", chrom, "with resolution of", binSize)
    minC = np.min(mat)
    #get the detected bins
    ss = {}
    for x, y in tqdm(mat):
        x = int((x - minC) / binSize)
        y = int((y - minC) / binSize)
        #if x == y:  #important to remove self-ligation
        #    continue
        if x not in ss:
            ss[x] = {}
        if y not in ss[x]:
            ss[x][y] = 0
        ss[x][y] += 1
    #only check non singleton bins
    keys = set()
    for x in list(ss.keys()):
        for y in list(ss[x].keys()):
            if ss[x][y] >= tol:
                keys.add((x, y))
    del ss
    #get the sub-sampling results
    rs = {}
    print("Get sub-sampling signals from", chrom)
    for ratio in tqdm(np.arange(0.05, 1, 0.05)):
        rs[ratio] = {}
        for rep in range(repeats):
            nmat = mat[np.random.choice(mat.shape[0], int(mat.shape[0] *
                                                          ratio)), :]
            ns = {}
            for x, y in nmat:
                x = int((x - minC) / binSize)
                y = int((y - minC) / binSize)
                #if x == y:  #important to remove self-ligation
                #    continue
                if x not in ns:
                    ns[x] = {}
                if y not in ns[x]:
                    ns[x][y] = 0
                ns[x][y] += 1
            c = 0
            #check how many bins still obtained
            for (x, y) in keys:
                if x in ns and y in ns[x] and ns[x][y] >= tol:
                    c += 1
            rs[ratio][rep] = c
    rs = pd.DataFrame(rs).T
    return mat.shape[0], len(keys), rs


def estSat(predir,
           fnOut,
           logger,
           bs=[10000, 5000],
           tol=5,
           cpu=1,
           cut=0,
           mcut=-1):
    """
    Estimation of sequencing saturation based on contact matrix signal ssampling capture.
    """
    metaf = predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    fs = []
    for key in meta["data"]["cis"]:
        fs.append(meta["data"]["cis"][key]["ixy"])
    samps = []
    for binSize in bs:
        ds = Parallel(n_jobs=cpu, backend='multiprocessing')(
            delayed(getSampBins)(f, binSize, cut=cut, mcut=mcut, tol=tol)
            for f in fs)
        totPETs = 0
        totBins = 0
        samp = []
        i = 0
        for d in ds:
            if d is not None:
                totPETs += d[0]
                totBins += d[1]
                if i == 0:
                    samp = d[2]
                    i = i + 1
                else:
                    samp = samp + d[2]
        samp = samp / totBins
        samps.append(samp)
    plotEstSat(bs, totPETs, samps, tol, fnOut)
