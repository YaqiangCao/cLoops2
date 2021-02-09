#!/usr/bin/env python
#--coding:utf-8 --
"""
filterWithBg.py
Get the filtered data against control data for 1D, then can use ixy2bdg.py to generate bedGraph file.
May filter too many
"""
__date__ = "2020-01-09"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os, argparse,json,random
from glob import glob
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#cLoops
from cLoops2.io import parseIxy
from cLoops2.settings import *


def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Filter the target data with control data (as substract).
        Example:
        filterWithBg.py -d CTCF_chic -bgd IgG_ChIC -o CTCF_ChIC_norm
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-d",
                        dest="dir",
                        required=True,
                        type=str,
                        help="Directory of target data from cLoops2 pre generated.")
    parser.add_argument(
        "-bgd",
        dest="bgdir",
        required=True,
        type=str,
        help="Directory of control data from cLoops2 pre generated.")
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    parser.add_argument('-p',
                        dest="cpu",
                        required=False,
                        default=1,
                        type=int,
                        help="Number of CPUs to run the job, default is 1.")
    parser.add_argument('-sf',
                        dest="sf",
                        required=False,
                        default=0.0,
                        type=float,
                        help="Scaling factor, by default will be estimated from data.")
    parser.add_argument('-bs',
                        dest="bs",
                        required=False,
                        default=150,
                        type=float,
                        help="Bin size to normalize the data, default is 150 for 1D track such as ChIC-seq/ChIP-seq/ATAC-seq.")
    op = parser.parse_args()
    return op


def getObs(xy,start,r):
    """
    Get the observed contact matrix in dict.
    """
    obsDict = {}
    for i,(x, y) in enumerate(xy):
        nx = int((x - start) / r)
        ny = int((y - start) / r)
        if nx not in obsDict:
            obsDict[nx] = {}
        if ny not in obsDict[nx]:
            obsDict[nx][ny] = {"counts":0,"rawId":[]}
        obsDict[nx][ny]["counts"] += 1
        obsDict[nx][ny]["rawId"].append( i )
    return obsDict


def filterWithBg(key,fixy,bgfixy,outdir,sf,r=150):
    obsKey,obsMat = parseIxy( fixy)
    obsXy = obsMat[:,1:]
    start = np.min(obsXy)
    obs = getObs( obsXy,start,r)

    expKey,expMat = parseIxy(bgfixy)
    expXy = expMat[:,1:]
    exp = getObs( expXy,start,r)
    
    rids = []
    #normalize the data
    for nx in exp.keys():
        if nx not in obs:
            continue
        for ny in exp[nx].keys():
            if ny not in obs[nx]:
                continue
            n = round(sf * exp[nx][ny]["counts"])
            if n > obs[nx][ny]["counts"]:
                del obs[nx][ny]
            else:
                rs = obs[nx][ny]["rawId"]
                random.shuffle(rs)
                rs = rs[n:]
                rids.extend(rs)
    print(obsMat.shape,len(rids))





def main():
    op = help()

    #prepare directory
    if not os.path.exists(op.dir):
        print("ERROR! %s not exists. return." % op.dir)
        return
    if not os.path.exists(op.bgdir):
        print("ERROR! %s not exists. return." % op.bgdir)
        return
    fdir = op.output
    if not os.path.exists(fdir):
        os.mkdir(fdir)
    elif len(os.listdir(fdir)) > 0:
        r = "ERROR! Working directory %s exists and not empty." % fdir
        print(r)
        return

    meta = json.loads(open(op.dir + "/petMeta.json").read())
    bgmeta = json.loads( open(op.bgdir + "/petMeta.json").read() )
    if op.sf == 0.0:
        sf = meta["Unique PETs"] / float(bgmeta["Unique PETs"])
    else:
        sf = op.sf

    keys = list(meta["data"]["cis"].keys())
    bgkeys = list(bgmeta["data"]["cis"].keys())
    for key in keys:
        if key not in bgkeys:
            print("ERROR! chrom %s in target data but not in control data. Please run cLoops2 pre with the same -c.")
            return
    
    #run the job
    Parallel(n_jobs=op.cpu)(
        delayed(filterWithBg)(
                           key,
                           meta["data"]["cis"][key]["ixy"],
                           bgmeta["data"]["cis"][key]["ixy"],
                           fdir,
                           sf,
                           op.bs,
                           ) for key in keys)

    
if __name__ == "__main__":
    main()
