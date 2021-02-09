#!/usr/bin/env python
#--coding:utf-8--
"""
cLoops2: dump
cLoops2 major file conversion functions
2020-03-04: to add new washU/UCSC support bigInteract track, according to https://genome.ucsc.edu/goldenPath/help/interact.html 
2020-06-25: to add dump to BEDPE
2020-07-01: refine ixy2bdg
2020-07-28: ixy2bed added
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import warnings
warnings.filterwarnings("ignore")
import os
import gzip 
import json
import random
from glob import glob

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm

#cLoops2
from cLoops2.io import parseIxy
from cLoops2.cmat import getObsMat, getExpMat
from cLoops2.utils import isTool, callSys
from cLoops2.settings import *

def ixy2bed(
            d, 
            fout, 
            logger,
            cut=0, 
            mcut=-1,
            ext=50,
    ):
    """
    Convert PETs to sorted BED file.
    @param d: str,cLoops2 pre data directory
    @param fout: str,prefix of output files
    @param logger: logger
    @param cut: int, > cut PETs kept
    @param mcut: int, <mcut PETs kept
    @param ext: int, extension from the PET center
    """
    if not os.path.exists(d):
        logger.error("%s not exists. return." % d)
        return
    fout = fout + "_reads.bed.gz"
    if os.path.isfile(fout):
        logger.error("Traget output file %s has been generated, return."%fout)
        return

    logger.info("Converting %s to BED file." % d)
    fs = glob(d + "/*.ixy")
    nfs = []
    for f in fs:
        chrom = f.split("/")[-1].split(".ixy")[0].split("-")
        if chrom[0] == chrom[1]:
            nfs.append(f)
    fs = nfs
    fs.sort()  #sort as bed files requries lexicograhic
    i = 0
    with gzip.open(fout, "wt") as f:
        for fin in fs:
            print("converting %s" % fin)
            key, mat = parseIxy(fin, cut=cut,mcut=mcut)
            s = list(mat[:,0])
            s.extend( list(mat[:,1]) )
            del mat
            s = np.array(list(set(s)))
            s.sort()
            for t in tqdm(s):
                line = [key[0], max([0, t - ext]), t + ext]
                f.write("\t".join(map(str, line)) + "\n")
    logger.info("Converting to BED file %s finished." % fout)


def ixy2bedpe(
            d, 
            fout, 
            logger,
            cut=0, 
            mcut=-1,
            ext=50,
    ):
    """
    Convert PETs to BEDPE file.
    @param d: str,cLoops2 pre data directory
    @param fout: str,prefix of output files
    @param logger: logger
    @param cut: int, > cut PETs kept
    @param mcut: int, <mcut PETs kept
    @param ext: int, extension from the PET center
    """
    if not os.path.exists(d):
        logger.error("%s not exists. return." % d)
        return
    fout = fout + "_PETs.bedpe.gz"
    if os.path.isfile(fout):
        logger.error("Traget output file %s has been generated, return."%fout)
        return

    logger.info("Converting %s to BEDPE file." % d)
    fs = glob(d + "/*.ixy")
    nfs = []
    for f in fs:
        chrom = f.split("/")[-1].split(".ixy")[0].split("-")
        if chrom[0] == chrom[1]:
            nfs.append(f)
    fs = nfs
    fs.sort()  #sort as bed files requries lexicograhic
    i = 0
    with gzip.open(fout, "wt") as f:
        for fin in fs:
            print("converting %s" % fin)
            key, mat = parseIxy(fin, cut=cut,mcut=mcut)
            for t in tqdm(mat):
                a = (key[0], max([0, t[0] - ext]), t[0] + ext)
                b = (key[1], max([0, t[1] - ext]), t[1] + ext)
                line = [
                    a[0], a[1], a[2],
                    b[0], b[1], b[2], 
                    i, ".", "+","-",
                ]
                f.write("\t".join(map(str, line)) + "\n")
                i += 1
    logger.info("Converting to BEDPE file %s finished." % fout)



def ixy2hic(
            d,
            fout,
            logger,
            org="hg38",
            resolution="1000,5000,10000,50000,100000,200000",
            cut=0,
            mcut=-1,
    ):
    """
    Convert reads level bedpe to HIC.
    Track format according to https://github.com/theaidenlab/juicer/wiki/Pre#file-format
    @param d: str,cLoops2 pre data directory
    @param fout: str,prefix of output files
    @param logger: logger
    @param org: str, species for the data
    @param resolution: str, series of resolutions
    @param cut: int, > cut PETs kept
    @param mcut: int, <mcut PETs kept
    """
    if not isTool("juicer_tools"):
        logger.error( "juicer_tools is not available in the executative enviroment! Please install and re-try.")
        return
    if not os.path.exists(d):
        logger.error("%s not exists,return." % d)
        return
    if os.path.isfile(fout+".hic"):
        logger.error("Traget output file %s.hic has been generated, return."%fout)
        return

    logger.info("Converting %s to .hic file which could be loaded in juicebox" % d)
    fs = glob(d + "/*.ixy")
    tmp = str(random.random())
    ss = {"+": 0, "-": 1}
    with open(tmp, "w") as fo:
        for fin in fs:
            print("converting %s" % fin)
            key, mat = parseIxy(fin, cut=cut,mcut=mcut)
            for t in tqdm(mat):
                line = [0, key[0], t[0], 0, 1, key[1], t[1], 1]
                fo.write("\t".join(map(str, line)) + "\n")
    # -n option from juicer_tools will cause trouble for trac-looping data
    #c1 = "juicer_tools pre -n -t ./ -r {resolution} {fin} {fout} {org}".format( 
    c1 = "juicer_tools pre -t ./ -r {resolution} {fin} {fout} {org}".format(
        resolution=resolution, fin=tmp, fout=fout+".hic", org=org)
    c2 = "rm %s" % tmp
    callSys([c1, c2])
    logger.info("Converting to juicer's hic file %s finished." % fout)



def ixy2washU(
            d, 
            fout, 
            logger,
            cut=0, 
            mcut=-1,
            ext=50,
    ):
    """
    Convert PETs to washU long range interactions. 
    Track format according to http://wiki.wubrowse.org/Long-range
    @param d: str,cLoops2 pre data directory
    @param fout: str,prefix of output files
    @param logger: logger
    @param cut: int, > cut PETs kept
    @param mcut: int, <mcut PETs kept
    @param ext: int, extension from the PET center
    """
    for t in ["bgzip", "tabix"]:
        if not isTool(t):
            logger.error(
                "%s is not available in the executative enviroment! Please install and re-try."
                % t)
            return
    if not os.path.exists(d):
        logger.error("%s not exists. return." % d)
        return
    fout = fout + "_PETs_washU.txt"
    if os.path.isfile(fout+".gz"):
        logger.error("Traget output file %s.gz has been generated, return."%fout)
        return

    logger.info("Converting %s to washU track." % d)
    fs = glob(d + "/*.ixy")
    nfs = []
    for f in fs:
        chrom = f.split("/")[-1].split(".ixy")[0].split("-")
        if chrom[0] == chrom[1]:
            nfs.append(f)
    fs = nfs
    fs.sort()  #sort as bed files requries lexicograhic
    i = 0
    with open(fout, "w") as f:
        for fin in fs:
            print("converting %s" % fin)
            key, mat = parseIxy(fin, cut=cut,mcut=mcut)
            #duplicate mat and swap column, for convient of sorting
            mat2 = np.copy(mat)
            mat2[:, [0, 1]] = mat2[:, [1, 0]]
            mat = np.concatenate((mat, mat2), axis=0)
            inds = np.argsort(mat[:, 0])
            mat = mat[inds, :]
            for t in tqdm(mat):
                a = (key[0], max([0, t[0] - ext]), t[0] + ext)
                b = (key[1], max([0, t[1] - ext]), t[1] + ext)
                line = [
                    a[0], a[1], a[2],
                    "%s:%s-%s,1" % (b[0], b[1], b[2]), i, "."
                ]
                f.write("\t".join(map(str, line)) + "\n")
                i += 1
    c1 = "bgzip %s" % fout
    c2 = "tabix -p bed %s.gz" % fout
    callSys([c1, c2])
    logger.info("Converting to washU long-range track %s finished." % fout)



def ixy2ucsc(
            d, 
            fout, 
            chromSizeF,
            logger,
            cut=0, 
            mcut=-1,
            ext=50,
    ):
    """
    Convert PETs to UCSC bigInteract track. 
    Track format according to https://genome.ucsc.edu/goldenPath/help/interact.html
    @param d: str,cLoops2 pre data directory
    @param fout: str,prefix of output files
    @param chromSizeF: str, file of chrom sizes, can be obtained through fetchChromSizes
    @param logger: logger
    @param cut: int, > cut PETs kept
    @param mcut: int, <mcut PETs kept
    @param ext: int, extension from the PET center

    """
    for t in ["bedToBigBed"]:
        if not isTool(t):
            logger.error(
                "%s is not available in the executative enviroment! Please install and re-try."
                % t)
            return
    if not os.path.exists(d):
        logger.error("%s not exists. return." % d)
        return
    if not os.path.exists(d):
        logger.error("%s not exists. return." % chromSizeF)
        return
    if os.path.isfile(fout+".bb"):
        logger.error("Traget output file %s.bb has been generated, return."%fout)
        return

    fildes="""
table interact
"Interaction between two regions"
    (
    string chrom;      "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records"
    uint chromStart;   "Start position of lower region. For interchromosomal, set to chromStart of this region"
    uint chromEnd;     "End position of upper region. For interchromosomal, set to chromEnd of this region"
    string name;       "Name of item, for display.  Usually 'sourceName/targetName' or empty"
    uint score;        "Score from 0-1000."
    double value;      "Strength of interaction or other data value. Typically basis for score"
    string exp;        "Experiment name (metadata for filtering). Use . if not applicable"
    string color;      "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4."
    string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
    uint sourceStart;  "Start position source/lower/this region"
    uint sourceEnd;    "End position in chromosome of source/lower/this region"
    string sourceName;  "Identifier of source/lower/this region"
    string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
    string targetChrom; "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
    uint targetStart;  "Start position in chromosome of target/upper/this region"
    uint targetEnd;    "End position in chromosome of target/upper/this region"
    string targetName; "Identifier of target/upper/this region"
    string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"
    )
    """
    #tmp file
    with open(fout+".tmp.as","w") as fo:
        fo.write(fildes)

    logger.info("Converting %s to UCSC track." % d)
    fs = glob(d + "/*.ixy")
    nfs = []
    for f in fs:
        chrom = f.split("/")[-1].split(".ixy")[0].split("-")
        if chrom[0] == chrom[1]:
            nfs.append(f)
    fs = nfs
    fs.sort()  #sort as bed files requries lexicograhic
    i = 0
    with open(fout+".tmp.bed", "w") as f:
        for fin in fs:
            print("converting %s" % fin)
            key, mat = parseIxy(fin, cut=cut,mcut=mcut)
            #duplicate mat and swap column, for convient of sorting
            mat2 = np.copy(mat)
            mat2[:, [0, 1]] = mat2[:, [1, 0]]
            mat = np.concatenate((mat, mat2), axis=0)
            inds = np.argsort(mat[:, 0])
            mat = mat[inds, :]
            for t in tqdm(mat):
                a = (key[0], max([0, t[0] - ext]), t[0] + ext)
                b = (key[1], max([0, t[1] - ext]), t[1] + ext)
                line = [
                    a[0], a[1], a[2],
                    ".",1,1,".",0,
                    a[0], a[1], a[2],".",".",
                    b[0], b[1], b[2], ".",".",
                ]
                f.write("\t".join(map(str, line)) + "\n")
                i += 1
    c1 = "bedToBigBed -tab -as=%s.tmp.as -type=bed5+13 %s.tmp.bed %s %s.bb"%( fout,fout,chromSizeF,fout )
    c2 = "rm %s.tmp.bed %s.tmp.as"%(fout,fout)
    callSys([c1,c2])
    logger.info("Converting to UCSC bigInteract track %s finished." % fout)



def addCov(cov, iv):
    """
    Add the coverage for a array. No value region is marked as False.
    """
    if len(cov) < iv[1]:
        cov.extend([False] * (iv[1] - len(cov) + 1))
    for i in range(iv[0], iv[1]):
        if cov[i] == False:
            cov[i] = 0
        cov[i] += 1
    return cov


def ixy2bdg(    
        d, 
        fout, 
        logger,
        cut=0, 
        mcut=-1,
        ext=50,
        pe=False,
    ):
    """
    Convert PETs to 1D bedGraph file with intrac-chromosomal PETs.

    @param d: str,cLoops2 pre data directory
    @param fout: str,prefix of output files
    @param logger: logger
    @param cut: int, > cut PETs kept
    @param mcut: int, <mcut PETs kept
    @param pe: bool, whether to treat paired-end tags as single end, set to True for ChIP-seq, ATAC-seq
    """
    if not os.path.exists(d):
        logger.error("%s not exists. return." % op.dir)
        return
    fout = fout + ".bdg"
    if os.path.isfile(fout):
        logger.error("Traget output file %s has been generated, return."%fout)
        return

    logger.info("Converting %s to bedGraph, normalized as RPM." % d)
    fs = glob(d + "/*.ixy")
    nfs = []
    for f in fs:
        chrom = f.split("/")[-1].split(".ixy")[0].split("-")
        #only cis PETs used
        if chrom[0] == chrom[1]:
            nfs.append(f)
    fs = nfs
    fs.sort()  #sort as lexicograhicly

    metaf = d + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    tot = meta["Unique PETs"] * 2
    
    with open(fout, "w") as f:
        for fin in fs:
            print("converting %s" % fin)
            key, mat = parseIxy(fin, cut=cut,mcut=mcut)
            #obtain coverage
            cov = []
            for t in tqdm(mat):
                if pe: #not interacting data, only normal paired-end, used for ChIC-seq, ATAC-seq.
                    a = (max([0, t[0] - ext]), t[1] + ext)
                    cov = addCov(cov, a)
                else:
                    a = (max([0, t[0] - ext]), t[0] + ext)
                    cov = addCov(cov, a)
                    b = (max([0, t[1] - ext]), t[1] + ext)
                    cov = addCov(cov, b)
            #write the bedGraph file use the step vector
            print("writting %s to bedGraph" % key[0])
            i = 0
            while i < len(cov) - 1:
                if cov[i] == False:  #find the non 0 start
                    i += 1
                    continue
                for j in range(i + 1, len(cov)):  #find the same value stop
                    if cov[j] != cov[i]:
                        break
                v = cov[i] / 1.0 / tot * 10**6
                v = "%.3f"%v
                line = [key[0], i, j - 1, v]
                f.write("\t".join(list(map(str, line))) + "\n")
                if j == len(cov) - 1:
                    break
                i = j
            del cov
    logger.info("Converting to bedGraph track %s finished." % fout)


def ixy2mat(
        d,
        fout,
        logger,
        chrom="",
        start=-1,
        end=-1,
        r=5000,
        cut=0,
        mcut=-1,
        log=False,
        method="obs",
        corr=False,
        norm=False,
):
    """
    Get the contact matrix.
    @param d: str,cLoops2 pre data directory
    @param fout: str,prefix of output files
    @param logger: logger
    @param chrom: str, such as "chr1-chr1"
    @param start: int, start location
    @param end: int, end location
    @param r: int, resolution bin size
    @param cut: int, > cut PETs kept
    @param mcut: int, <mcut PETs kept
    @param log: bool, whether do log transformation
    @param method: choice, available are obs, obs/exp
    @param corr: bool, whehter to get the correlation matrix
    @param norm: bool, whether to normalize the matrix with z-score
    """
    if start != -1 and end != -1 and end < start:
        logger.error("End %s is smaller than %s start." % (end, start))
        return
    f = os.path.join(d,chrom+".ixy")
    if not os.path.isfile(f):
        logger.error("%s not exists, please check the input -mat_chrom"%f)
        return

    chrom, xy = parseIxy(f, cut=cut,mcut=mcut)
    if start == -1:
        start = np.min(xy)
    if end == -1:
        end = np.max(xy)
    mat = getObsMat(xy, start, end, r)
    bgmat = None
    if method == "obs/exp":
        bgmat = getExpMat(xy, mat.shape, start, end, r)
    if log:
        if bgmat is None:
            mat = np.log10(mat + 1)
        else:
            mat = np.log10(mat + 1) - np.log10(bgmat + 1)
    else:
        if bgmat is not None:
            mat = (mat + 1) / (bgmat + 1)
    if corr:
        mat = np.corrcoef(mat)
        mat = np.nan_to_num(mat)
    if norm:
        m = np.mean(mat)
        s = np.std(mat)
        mat = (mat - m) / s
    rs = []
    bs = int((end - start) / r) + 1
    for i in range(bs):
        nr = (chrom[0], start + i * r, start + i * r + r)
        rs.append("|".join(list(map(str, nr))))
    mat = pd.DataFrame(mat, index=rs, columns=rs)
    mat.to_csv(fout + "_cmat.txt", sep="\t", index_label="pos")
    logger.info("Converting to contact matrix txt %s_cmat.txt finished." % fout)

