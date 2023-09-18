#!/home/caoy7/anaconda3/envs/astroBoy/bin/python
#--coding:utf-8 --
"""
getFRiF.py
Get the fraction of PETs in target features. 
"""
__date__ = "2020-01-09"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os, argparse, json
from glob import glob
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

#cLoops
from cLoops2.ds import XY, Loop
from cLoops2.io import parseIxy, parseTxt2Loops
from cLoops2.settings import *


def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Get the fractrion of reads in target features such as domains 
        and peaks annotated with .bed file or domains/stripes/loops with 
        .txt file such as the _loop.txt file.

        getFRiF.py -d GM12878_Trac -b tad.bed 
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-d",
                        dest="predir",
                        required=True,
                        type=str,
                        help="Directory for cLoops2 pre generated.")
    parser.add_argument(
        "-b",
        dest="bed",
        required=False,
        default="",
        type=str,
        help="The .bed annotated the features such peaks/domains.")
    parser.add_argument("-l",
                        dest="floop",
                        required=False,
                        type=str,
                        default="",
                        help="The _loop.txt .file generated by cLoops2.")
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    parser.add_argument(
        "-single",
        dest="single",
        required=False,
        default=False,
        action="store_true",
        help=
        "Whether to treat paired-end as single when check the location in -bed features, default is False, which means the two ends has to be in the same target regions, then it counts one. Set this when -bed features are peaks."
    )

    parser.add_argument(
        "-cut",
        dest="cut",
        type=int,
        default=0,
        help=
        "Distance cutoff for PETs to filter, default is 0. Can be set as the estimated self-ligation distance cutoff."
    )
    parser.add_argument('-p',
                        dest="cpu",
                        required=False,
                        default=1,
                        type=int,
                        help="Number of CPUs to run the job, default is 1.")
    op = parser.parse_args()
    return op


def get1DFRiF(key, fixy, peaks, cut=0, single=False):
    """
    Get reads ratio within 1D peaks/domains.
    """
    key2, mat = parseIxy(fixy, cut=cut)
    if mat.shape[0] == 0:
        print(
            "No PETs found in %s maybe due to distance cutoff for PET > %s." %
            (fixy, cut))
    xy = XY(mat[:, 0], mat[:, 1])
    if single:
        print("Getting the single end tags target %s 1d features from %s." %
              (len(peaks), key))
    else:
        print("Getting the paired end tags target %s 1d features from %s." %
              (len(peaks), key))
    ns = 0
    for peak in tqdm(peaks):
        if single:
            ns += len(xy.queryPeak(peak[1], peak[2]))
        else:
            ns += len(xy.queryPeakBoth(peak[1], peaks[2]))
    if single:
        print("%s: total %s single end tags, %s (%.3f) in target." %
              (key, mat.shape[0] * 2, ns, float(ns) / mat.shape[0] / 2))
    else:
        print("%s: total %s paired end tags, %s (%.3f) in target." %
              (key, mat.shape[0], ns, float(ns) / mat.shape[0]))
    return ns


def get2DFRiF(key, fixy, loops, cut=0):
    """
    Get reads ratio within 2D loops/domains.
    """
    key2, mat = parseIxy(fixy, cut=cut)
    if mat.shape[0] == 0:
        print(
            "No PETs found in %s maybe due to distance cutoff for PET > %s." %
            (fixy, cut))
    xy = XY(mat[:, 0], mat[:, 1])
    ns = 0
    for loop in tqdm(loops):
        ns += len(
            xy.queryLoop(loop.x_start, loop.x_end, loop.y_start,
                         loop.y_end)[-1])
    print("%s: total %s paired end tags, %s (%.3f) in target." %
          (key, mat.shape[0], ns, float(ns) / mat.shape[0]))
    return ns


def main():
    op = help()
    if op.bed == "" and op.floop == "":
        print("No 1D or 2D features assigned as input, return.")
        return
    metaf = op.predir + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    keys = list(meta["data"]["cis"].keys())
    keys.extend(list(meta["data"]["trans"].keys()))
    if op.bed != "" and os.path.exists(op.bed):
        peaks = {}
        for line in open(op.bed):
            line = line.split("\n")[0].split("\t")
            if len(line) < 3:
                continue
            key = line[0] + "-" + line[0]
            if key not in peaks:
                peaks[key] = []
            line[1] = int(line[1])
            line[2] = int(line[2])
            peaks[key].append(line)
        ds = Parallel(n_jobs=op.cpu,backend="multiprocessing")(delayed(get1DFRiF)(
            key,
            meta["data"]["cis"][key]["ixy"],
            peaks[key],
            cut=op.cut,
            single=op.single,
        ) for key in keys if key in peaks)
        #ds = [d[0] for d in ds]
        s = np.sum(ds)
        t = meta["Unique PETs"]
        if op.single:
            t = t * 2
        ds = {
            "total": t,
            "inTraget": s,
            "ratio": float(s) / t,
            "single": op.single
        }
        ds = pd.Series(ds)
        ds.to_csv(op.output + "_1d_FRiF.txt", sep="\t")

    if op.floop != "" and os.path.exists(op.floop):
        loops = parseTxt2Loops(op.floop)
        ds = Parallel(n_jobs=op.cpu,backend="multiprocessing")(delayed(get2DFRiF)(
            key,
            meta["data"]["cis"][key]["ixy"],
            loops[key],
            cut=op.cut,
        ) for key in keys if key in loops)
        #ds = [d[0] for d in ds]
        s = np.sum(ds)
        t = meta["Unique PETs"]
        ds = {
            "total": t,
            "inTraget": s,
            "ratio": float(s) / t,
        }
        ds = pd.Series(ds)
        ds.to_csv(op.output + "_2d_FRiF.txt", sep="\t")


if __name__ == "__main__":
    main()
