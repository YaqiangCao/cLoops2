#!/usr/bin/env python3
#--coding:utf-8 --

#sys
from math import ceil

#3rd
import numpy as np
import pandas as pd
from scipy import linalg
from scipy.stats import chi2


def lowess(x, y, f=2. / 3., iter=3):
    """lowess(x, y, f=2./3., iter=3) -> yest

    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.

    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.
    """
    n = len(x)
    r = int(ceil(f * n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    w = (1 - w**3)**3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
            A = np.array([[np.sum(weights),
                           np.sum(weights * x)],
                          [np.sum(weights * x),
                           np.sum(weights * x * x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]

        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta**2)**2
    return yest


def getBonPvalues(ps):
    """
    Return the Bonferroni corrected p-values.
    """
    ps = np.array(ps)
    ps = ps * len(ps)
    ps[ps > 1.0] = 1.0
    return ps


def mahalanobis(mat):
    """
    Caculate the mahalanobis distance.

    according to: https://www.statology.org/mahalanobis-distance-python/

    @param mat: np.array, row is each item and column is each variable 

    @return dis: np.array,Mahalanobis distance
    @return ps: np.array, chi-square test p-values
    """
    #covariance matrix 
    cov = np.cov(mat, rowvar=False)
    #inverse covariance matrix 
    invCov = np.linalg.inv(cov)
    #center 
    center = np.mean(mat, axis=0)
    #mahalanobis distance
    mu = mat - center
    dis = np.dot( np.dot(mu,invCov), mu.T ).diagonal()
    #Chi-square test p-values for detecting outliers 
    ps = 1 - chi2.cdf(dis, mat.shape[1]-1)
    return dis, ps


def twoPassesMDTest(data, pcut=0.01):
    """
    Perform MD test with two passes. First pass use all data, second pass using cov and center from data without outliers. 
    @param data: pd.DataFrame, row is item and column is sample.
    @param pcut: Chi-Square p-value cutoffs
    """
    #first pass test
    #mahalanobis distance and pvaues
    dis, ps = mahalanobis(data.values)
    dis = pd.Series(dis, index=data.index)
    ps = pd.Series(ps, index=data.index)
    inds = ps[ps < pcut].index

    #second pass test, only using data without outliers to caculate cov matrix and center
    ndata = data.drop(inds)
    #covariance matrix
    cov = np.cov(ndata.values, rowvar=False)
    #inverse covariance matrix
    invCov = np.linalg.inv(cov)
    #center
    center = np.mean(ndata.values, axis=0)

    #mahalanobis distance for all data
    mu = data.values - center
    dis = np.dot(np.dot(mu, invCov), mu.T).diagonal()
    #Chi-square test p-values for detecting outliers from all data
    ps = 1 - chi2.cdf(dis, data.shape[1] - 1)

    dis = pd.Series(dis, index=data.index)
    ps = pd.Series(ps, index=data.index)
    return dis, ps


