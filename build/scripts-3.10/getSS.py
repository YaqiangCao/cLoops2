#!/home/caoy7/anaconda3/envs/astroBoy/bin/python
#--coding:utf-8 --
"""
getSS.py
cLoops2 getSS.py caculation the segregation score for domain-centric analysis.

"""
__date__ = "2022-10-17"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os, argparse
from glob import glob
from collections import Counter
from datetime import datetime
from argparse import RawTextHelpFormatter

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



def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Caculate the segregation score for a specific region. 
        The output .bdg is the bedGraph result for the regions with segregation score.

        Example:
        getSS.py -f GM12878_Trac/chr21-chr21.ixy -o GM12878_Trac_chr21
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-f",
        dest="fixy",
        required=True,
        type=str,
        help=
        "Input .ixy file generated by cLoops2 to caculate segregation score.")
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    parser.add_argument(
        "-bs",
        dest="binSize",
        required=False,
        default=1000,
        type=int,
        help="Bin size (bp) to generate the contact matrix, default is 1000 bp."
    )
    parser.add_argument(
        "-ws",
        dest="winSize",
        required=False,
        default=50000,
        type=int,
        help=
        "The half of the sliding window size used to caculate local correlation, default is 50000 (50kb)."
    )
    parser.add_argument(
        "-cut",
        dest="cut",
        type=int,
        default=0,
        help="PETs with distance > cut will be kept, default is 0.")
    parser.add_argument(
        "-mcut",
        dest="mcut",
        type=int,
        default=-1,
        help="PETs with distance < mcut will be kept, default is -1 no limit.")
    op = parser.parse_args()
    return op



def calcSS(f, fout, bs=20000, winSize=500000, cut=0,mcut=-1):
    """
    Calculation of segregation score, output as .bedGraph file.
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
        nmat[nmat < 0] = 0
        s = nmat.mean()
        ss.append(s)
        r = [key[0], x, x + bs]
        rs.append(r)
    ss = np.array(ss)
    ss = (ss - np.mean(ss))/np.std(ss)
    for i, r in enumerate(rs):
        r.append(ss[i])
    with open(fout +"_SS.bdg", "w") as fo:
        for r in rs:
            fo.write("\t".join(list(map(str, r))) + "\n")




def main():
    op = help()
    calcSS(op.fixy,
           op.output,
           bs=op.binSize,
           winSize=op.winSize,
           cut=op.cut,
           mcut=op.mcut,
           )


if __name__ == "__main__":
    main()
