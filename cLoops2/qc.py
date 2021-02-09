#!/usr/bin/env python
#--coding:utf-8--

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
import gc
import sys
import time
import gzip
from datetime import datetime

#3rd library
from joblib import Parallel, delayed

#cLoops2
from cLoops2.ds import PET
from cLoops2.utils import cFlush


class bedpeSTAT:
    """
    basic data structure for a bedpe file stat
    """
    sample = None
    total = 0  #total PETs
    cis = 0  #cis PETs
    close = 0  #distance < 1k
    mid = 0  #distance > 1k < 10k
    far = 0  #distance > 10k
    fr = 0  #one mapped to postive and another mapped to negative strand
    ff = 0
    rr = 0
    uniques = 0
    trans = 0
    redudancy = 0
    meandis = 0
    #chrMRatio = 0


def evaBedpe(f):
    """
    Qaulity control for the bedpe file, the bedpe file is coverted from .bam file, contains all information.
    """
    stat = bedpeSTAT()
    reds = set()
    if f.endswith(".gz"):
        of = gzip.open(f, "rt")  #python 3
    else:
        of = open(f)
    print("Start parsing %s" % f)
    for i, line in enumerate(of):
        if i % 100000 == 0:
            cFlush("%s PETs processed from %s" % (i, f))
        line = line.split("\n")[0].split("\t")
        if len(line) < 6:
            continue
        try:
            pet = PET(line)
        except:
            continue
        stat.total += 1
        t = (pet.chromA, pet.chromB, pet.startA, pet.endA, pet.startB,
             pet.endB, pet.strandA, pet.strandB)
        t = hash(t)
        if t not in reds:
            reds.add(t)
        else:
            continue
        if pet.cis != True:
            stat.trans += 1
            continue
        #print(t)
        #only counting intra-chromosomal interaction PETs
        stat.cis += 1
        if pet.distance <= 1000:
            stat.close += 1
        if 1000 < pet.distance <= 10000:
            stat.mid += 1
        if 10000 < pet.distance:
            stat.far += 1
        if pet.strandA == "+" and pet.strandB == "-":
            stat.fr += 1
        if pet.strandA == "-" and pet.strandB == "+":
            stat.fr += 1
        if pet.strandA == "+" and pet.strandB == "+":
            stat.ff += 1
        if pet.strandA == "-" and pet.strandB == "-":
            stat.rr += 1
        stat.meandis += pet.distance
        #if pet.chromA == "chrM":
        #    chrMRatio += 1
    print()
    stat.sample = f.split("/")[-1].split(".bedpe")[0]
    stat.uniques = 0
    stat.uniques = len(reds)
    if stat.total > 0:
        stat.redudancy = 1 - stat.uniques / 1.0 / stat.total
    else:
        stat.redundancy = 0
    stat.meandis = stat.meandis / 1.0 / stat.cis
    print("Analysis of %s finished." % f)
    del reds
    gc.collect()
    return stat


def qcBedpes(fs, fout, cpu=1):
    """
    Quality control.
    """
    #run the job
    data = Parallel(n_jobs=min(len(fs), cpu),backend="multiprocessing")(delayed(evaBedpe)(f) for f in fs)

    with open(fout, "w") as fo:
        header = [
            "Sample",
            "TotalPETs",
            "UniquePETs",
            "Redundancy",
            "IntraChromosomalPETs(cis)",
            "cisRatio",
            "InterChromosomalPETs(trans)",
            "transRatio",
            "meanDistance",
            "closePETs(distance<=1kb)",
            "closeRatio",
            "middlePETs(1kb<distance<=10kb)",
            "middleRatio",
            "distalPETs(distance>10kb)",
            "distalRatio",
        ]
        fo.write("\t".join(header) + "\n")
        for bs in data:
            line = [
                bs.sample, bs.total, bs.uniques, bs.redudancy, bs.cis,
                bs.cis / 1.0 / bs.uniques, bs.trans,
                bs.trans / 1.0 / bs.uniques, bs.meandis, bs.close,
                bs.close / 1.0 / bs.cis, bs.mid, bs.mid / 1.0 / bs.cis, bs.far,
                bs.far / 1.0 / bs.cis
            ]
            line = list(map(str, line))
            fo.write("\t".join(line) + "\n")
