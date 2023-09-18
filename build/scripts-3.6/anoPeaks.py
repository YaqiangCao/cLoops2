#!/home/caoy7/anaconda2/envs/cLoops2/bin/python
#--coding:utf-8 --
"""
anoPeaks.py
cLoops2 anoPeaks.py annotate peaks genomic locations as promoters or enhancers.
"""

__date__ = "2023-08-10"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os
import argparse
from argparse import RawTextHelpFormatter

#3rd library
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.spatial import KDTree

#cLoops2
from cLoops2.ano import readGenes, findOverlapOrNearest
from cLoops2.ds import Peak


def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Annotate peaks as promoters or enhancers according to gene annotations. 

        Example:
        anoPeaks.py -f H3K27ac_peaks.bed -gtf hg38.gtf -o H3K27ac_peaks 
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-f",
        dest="peakf",
        required=True,
        type=str,
        help=
        "Input .bed file as peaks. Peak id will be renamed as chrom|start|end")
    parser.add_argument("-gtf",
                        dest="gtf",
                        default="",
                        required=False,
                        type=str,
                        help="GTF file annotation for genes.")
    parser.add_argument(
        "-tid",
        dest="tid",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to use transcript id instead of gene id for annotation. Default\n"\
        "is not."
    )
    parser.add_argument(
        "-pdis",
        dest="pdis",
        default=2000,
        required=False,
        type=int,
        help=
        "Distance limitation for anchor to nearest gene/transcript TSS to define\n"\
        "as promoter. Default is 2000 bp."
    )
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    op = parser.parse_args()
    return op


def parseBed2Peaks(fbed):
    """
    """
    peaks = {}
    for line in open(fbed):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        chrom = line[0]
        if chrom not in peaks:
            peaks[chrom] = []
        peak = Peak()
        peak.id = "|".join(line[:3])
        peak.chrom = chrom
        peak.start = int(line[1])
        peak.end = int(line[2])
        peak.length = peak.end - peak.start
        peaks[chrom].append(peak)
    return peaks


def anoPeaks(
        peakf,
        fout,
        gtf,
        tid=False,
        pdis=2000,
):
    """
    Annotate peaks.
    @param peakf: str, name of peaks file,  .bed file
    @param fout: str, output prefix
    @param gtf: str, GTF file name 
    @param tid: bool, if set true, use transcript id for alternative TSS
    @param pdis: <=distance nearest TSS to define as promoter
    """
    if not os.path.isfile(peakf):
        print("Input %s not exists, return." % peakf)
        return
    elif not os.path.isfile(gtf):
        print("Input %s not exists, return." % gtf)
        return
    else:
        #gene annotions, {chrom:{tss:g}}, tss is int
        genes = readGenes(gtf, tid=tid)
        peaks = parseBed2Peaks(peakf)
        #store results
        rs = {}
        #find nearest TSS
        print("annotating peaks" )
        for chrom in tqdm(peaks.keys()):
            if chrom not in genes:
                continue
            gs = genes[chrom]
            ts = np.array([[tss] for tss in gs.keys()])
            cov = {}
            for tss, g in gs.items():
                cov[tss] = g
            tree = KDTree(ts)
            for peak in peaks[chrom]:
                xgs, xds = findOverlapOrNearest(gs, ts, tree, peak.start,
                                                peak.end)
                if len(xgs) > 1:
                    xt = "Promoter"
                    xd = 0
                else:
                    xd = xds[0]
                    if abs(xd) <= pdis:
                        xt = "Promoter"
                    else:
                        xt = "Enhancer"
                rs[peak.id] = {
                    "1_chrom":
                    peak.chrom,
                    "2_start":
                    peak.start,
                    "3_end":
                    peak.end,
                    "4_type":
                    xt,
                    "5_nearestDistanceToTSS":
                    xd,
                    "6_nearestTargetTSS":
                    ",".join([
                        xg.chrom + ":" + str(xg.start) + "-" + str(xg.end) +
                        "|" + xg.strand + "|" + xg.name for xg in xgs
                    ]),
                }
        rs = pd.DataFrame(rs).T
        rs.to_csv(fout + "_anoPeaks.txt", sep="\t", index_label="peakId")


def main():
    op = help()
    anoPeaks(
        op.peakf,
        op.output,
        op.gtf,
        tid=op.tid,
        pdis=op.pdis,
    )


if __name__ == "__main__":
    main()
