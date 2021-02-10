#!python
#--coding:utf-8--
"""
tracPre.py
Pre-processing code for Hi-Trac data, implemented with cLoops2, from fastq to bedpe files and qc report.
2020-02-27: finished and well tested.
2020-06-30: add linker filter, new stat, and changing mapping to end-to-end
"""

__author__ = "CAO Yaqiang"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
import time
import gzip
import argparse
import subprocess
from glob import glob
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import pandas as pd
from joblib import Parallel, delayed
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#cLoops2
from cLoops2.utils import getLogger, callSys, isTool

#global settings
#logger
date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
logger = getLogger(fn=os.getcwd() + "/" + date.strip() + "_" +
                   os.path.basename(__file__) + ".log")


def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Preprocess the raw reads of FASTQ files of Trac-looping to reference
        geneome with bowtie2 and obtain the unqiue PETs with quality control
        results.
        Fastqs files should be named with suffix pattern as 
        _R1.fastq.gz, _R2.fastq.gz.

        Example:
        tracPre.py -fqd ../1.fq -o ./ -ref ../bowtie2/hg38 -n 10 -p 5 -mapq 10
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-fqd",
        dest="fqd",
        required=True,
        type=str,
        help="The directory for raw .fastq.gz files, for example ../1.fastq/ ")
    parser.add_argument(
        "-o",
        dest="output",
        required=False,
        type=str,
        default="./",
        help=
        "Output directory, default is ./, if directory not exists, create one."
    )
    parser.add_argument(
        "-ref",
        dest="ref",
        required=True,
        type=str,
        help=
        "Bowtie2 reference index prefix, such as ./ref/hg38, generated from\n"\
        "bowtie2-build hg38.fa hg38."
    )
    parser.add_argument(
        "-n",
        dest="number",
        required=False,
        type=int,
        default=1,
        help="How many Bowtie2 to run at the same time, default is 1. ")
    parser.add_argument(
        "-p",
        dest="cpu",
        required=False,
        type=int,
        default=5,
        help="How many cpus used by each Bowtie2 or following processing,\n"\
        "default is 5. "
    )
    parser.add_argument("-mapq",
                        dest="mapq",
                        required=False,
                        default=10,
                        type=int,
                        help="MAPQ cutoffs for filtering PETs, default is 10.")
    op = parser.parse_args()
    return op


def preFqs(fastqRoot):
    """
    If the fastq files are well prepared, suitable. 
    """
    fastqs = glob(fastqRoot + "/*.fastq.gz")
    data = {}
    for fq in fastqs:
        s = os.path.split(fq)[1]
        s = s.replace(".fastq.gz", "")
        if s.endswith("_R1"):
            sample = s.replace("_R1", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][0] = fq
        if s.endswith("_R2"):
            sample = s.replace("_R2", "")
            if sample not in data:
                data[sample] = [0, 0]
            data[sample][1] = fq
    for key, fqs in data.items():
        if len(fqs) != 2:
            logger.error(
                "for %s there is not paired fastq files, only %s found" %
                (key, ",".join(fqs)))
            del data[key]
    return data


def findLinker(seq, linker):
    """
    Match the linker in the read sequence.
    """
    pos = -1
    for i in range(len(seq) - 9):
        seed = seq[i:i + 9]
        if linker.startswith(seed):
            pos = i
            break
    return pos


def checkStarts(seq):
    """
    Check the starts 
    """
    flag = False
    ss = ["CATG", "AATT", "NATG", "NATT"]
    for s in ss:
        if seq.startswith(s):
            flag = True
            break
    return flag


def cutLinker(fq1, fq2, pre, rlen=10, linker="CTGTCTCTTATACACATCT"):
    """
    Cut linkers and filter too short reads
    """
    sample = pre.split("/")[-1]
    nf1 = pre + "_R1.fastq.gz"
    nf2 = pre + "_R2.fastq.gz"
    if os.path.isfile(nf1) and os.path.isfile(nf2):
        print("%s has been generated, return" % pre)
        return None
    fouts = {
        "fo_r1": gzip.open(nf1, "wt"),
        "fo_r2": gzip.open(nf2, "wt"),
    }
    #processing pairing fastqs
    with gzip.open(fq1, "rt") as f1, gzip.open(fq2, "rt") as f2:
        i = 0
        j = 0
        for r1, r2 in zip(FastqGeneralIterator(f1), FastqGeneralIterator(f2)):
            r1, r2 = list(r1), list(r2)
            i += 1
            if i % 100000 == 0:
                print("%s reads processed for %s" % (i, pre))
            #check the starts
            """
            if not (checkStarts(r1[1]) and checkStarts(r2[1])):
                continue
            if r1[1][0] == "N":
                r1[1] = r1[1][1:]
                r1[2] = r1[2][1:]
            if r2[1][0] == "N":
                r2[1] = r2[1][1:]
                r2[2] = r2[2][1:]
            """
            #check the linker
            r1pos = findLinker(r1[1], linker)
            r2pos = findLinker(r2[1], linker)
            #trim reads
            if r1pos != -1:
                r1[1] = r1[1][:r1pos]
                r1[2] = r1[2][:r1pos]
            if r2pos != -1:
                r2[1] = r2[1][:r2pos]
                r2[2] = r2[2][:r2pos]
            rid = "_".join(list(map(str, [i, r1pos, r2pos])))
            r1[0] = rid
            r2[0] = rid
            if len(r1[1]) >= rlen and len(r2[1]) >= rlen:
                j += 1
                fouts["fo_r1"].write("@%s\n%s\n+\n%s\n" %
                                     (r1[0], r1[1], r1[2]))
                fouts["fo_r2"].write("@%s\n%s\n+\n%s\n" %
                                     (r2[0], r2[1], r2[2]))
    return sample, i, j, nf1, nf2


def tracMapping(sample, fqs, ref, outdir, cpus=25):
    """
    Mapping settings for Trac-looping data.
    """
    logger.info("Start mapping %s.\n" % sample)
    od = os.path.join(outdir, sample)
    if not os.path.exists(od):
        os.makedirs(od, exist_ok=True)
    sam = od + "/" + sample + ".sam"
    bam = od + "/" + sample + ".bam"
    if os.path.isfile(sam):
        logger.error("%s:%s exists, return." % (sample, sam))
        return None
    if os.path.isfile(bam):
        logger.error("%s:%s exists, return." % (sample, bam))
        return None
    doBowtie = "bowtie2 -p {cpus} -q --end-to-end --very-sensitive -x {ref} -1 {fq1} -2 {fq2} -S {sam}".format(
        cpus=cpus, ref=ref, fq1=fqs[0], fq2=fqs[1], sam=sam)
    logger.info(doBowtie)
    stat, output = subprocess.getstatusoutput(doBowtie)
    #trim with "Warning"
    output = output.split("\n")
    output = [t for t in output if not t.startswith("Warning")]
    output = "\n".join(output)
    logger.info("FLAG_A:" + sample + "\n" + output + "\nFLAG_A\n")
    return sample, sam


def getUniqueBedpe(f, fout):
    """
    Get unique bedpe. Read id indicate the linker location.
    """
    if os.path.isfile(fout):
        return
    print("Getting unique PETs from %s to %s" % (f, fout))
    redus = set()
    with open(fout, "w") as fo:
        for i, line in enumerate(open(f)):
            line = line.split("\n")[0].split("\t")
            if len(line) < 6:
                continue
            rid = list(map(int, line[6].split("_")))
            #for cis short reads, requiring the linkers
            if line[0] == line[3]:
                dis = abs((int(line[1]) + int(line[2])) / 2 -
                          (int(line[4]) + int(line[5])) / 2)
                if dis < 1000 and rid[1] + rid[2] == -2:
                    continue
            #for trans reads, requiring the linkers
            if line[0] != line[3]:
                if rid[1] + rid[2] == -2:
                    continue
            #remove redudant PETs
            r = hash(tuple(line[:6]))
            if r in redus:
                continue
            else:
                redus.add(r)
            fo.write("\t".join(line) + "\n")


def sam2bamBedpe(sample, sam, mapq=10):
    """
    SAM to BAM and bedpe file 
    """
    n = os.path.splitext(sam)[0]
    bam = n + ".bam"
    bedpeAll = n + "_all.bedpe"
    bedpeUni = n + "_unique.bedpe"
    #sam to bam, filtering mapq
    samview = "samtools view -b -F 4 -@ 2 -q {mapq} -o {bam} {sam}".format(
        mapq=mapq, bam=bam, sam=sam)
    #sort by read name
    samsort = "samtools sort -n -@ 2 {bam} -T {pre} -o {bam}".format(
        bam=bam, pre=bam.replace(".bam", ""))
    rmsam = "rm %s" % (sam)
    cmds = [samview, samsort, rmsam]
    callSys(cmds, logger)
    bam2bedpe = "bamToBed -bedpe -i {bam} > {bedpe}".format(bam=bam,
                                                            bedpe=bedpeAll)
    logger.info(bam2bedpe)
    stat, output = subprocess.getstatusoutput(bam2bedpe)
    getUniqueBedpe(bedpeAll, bedpeUni)
    cmd = "gzip %s %s" % (bedpeAll, bedpeUni)
    callSys([cmd], logger)
    return sample, bedpeAll + ".gz", bedpeUni + ".gz"


def sParseBowtie(lines):
    """
    Parse Bowtie2 log file, to obtain mapping stastics.
    """
    d, s = None, None
    lines = lines.split("\n")
    s = lines[0]
    totalReads = int(lines[1].split(";")[0].split()[0])
    d1 = lines[4].strip().split()
    conUniqueMappedReads = int(d1[0])
    d2 = lines[8].strip().split()
    unconUniqueMappedReads = int(d2[0])
    #mapRatio = float(lines[15].split("%")[0])
    mapRatio = float(lines[-2].split("%")[0])
    d = {
        "TotalRawReads": totalReads,
        #"ConcordantlyUniqueMapReads": conUniqueMappedReads,
        #"DisconcordantlyUniqueMapReads": unconUniqueMappedReads,
        "MappingRatio(%s)": mapRatio
        #"MultipleMapReads": multipleMappedReads,
        #"MultipleMapRatio": multipleMappedRatio,
    }
    return d, s


def parseBowtielog(logs=None):
    if logs == None:
        logs = glob("*.log")
    data = {}
    for log in logs:
        lines = open(log).read().split("FLAG_A\n")
        lines = [line for line in lines if "FLAG_A" in line]
        for line in lines:
            t = line.split("FLAG_A:")[1]
            d, s = sParseBowtie(t)
            data[s] = d
    data = pd.DataFrame(data).T
    return data


def main():
    """
    Batch converting from bam to bedpe.
    """
    #prepare everything
    op = help()
    for t in ["bowtie2", "samtools", "bamToBed"]:
        if not isTool(t):
            logger.error("%s not exits! Please install through conda." % t)
            return
    if not os.path.exists(op.fqd):
        logger.error("Input %s not exists! Return." % op.fqd)
        return
    if len(glob(op.ref + "*.bt2")) == 0:
        logger.error("Bowtie2 reference not exists for prefix of %s! Return." %
                     op.ref)
        return
    if not os.path.exists(op.output):
        os.makedirs(op.output, exist_ok=True)
    else:
        fs = glob(os.path.join(op.output, "*"))
        if len(fs) > 0:
            logger.info(
                "Target output directory %s is not empty, may over-write some files."
                % op.output)
    data = preFqs(op.fqd)
    if len(data) == 0:
        logger.error(
            "No matched _R1.fastq.gz and _R2.fastq.gz in %s. Return." %
            (op.fqd))
        return
    #prepare output dir
    dirs = {}
    for sample in data.keys():
        od = os.path.join(op.output, sample)
        dirs[sample] = od
        if not os.path.exists(od):
            os.makedirs(od, exist_ok=True)

    #step 1, filter linkers
    logger.info("Step1: Trim linkers and remove short sequences.")
    ds = Parallel(n_jobs=op.number)(
        delayed(cutLinker)(fqs[0], fqs[1], os.path.join(dirs[sample], sample))
        for sample, fqs in data.items())
    data = {}
    for d in ds:
        if d is not None:
            data[d[0]] = {
                "totalRaw": d[1],
                "filterLinkers": d[2],
                "f1": d[3],
                "f2": d[4],
            }

    #step2, mapping
    logger.info("Step2: Map processed reads to genome.")
    ref = op.ref
    ds = Parallel(n_jobs=op.number, backend="multiprocessing")(
        delayed(tracMapping)(
            sample, [vs["f1"], vs["f2"]], ref, op.output, cpus=op.cpu)
        for sample, vs in data.items())
    for d in ds:
        if d is not None:
            data[d[0]]["sam"] = d[1]

    #step3, convert to bam and bedpe files
    #sam to bam and bedpe
    logger.info("Step3: File type conversion. ")
    cpus = op.number * op.cpu
    ncpus = int(min(len(data), cpus / 2))
    ds = Parallel(n_jobs=ncpus, backend="multiprocessing")(
        delayed(sam2bamBedpe)(sample, vs["sam"], op.mapq)
        for sample, vs in data.items())

    allBedpes = []
    uniBedpes = []
    for d in ds:
        if d is not None:
            data[d[0]]["allBedpe"] = d[1]
            data[d[0]]["uniNonbgBedpe"] = d[2]
            allBedpes.append(d[1])
            uniBedpes.append(d[2])
    data = pd.DataFrame(data).T

    #step 4, all PETs cLoops2 qc
    logger.info("Step4: All mapped PETs QC. ")
    cmd = "cLoops2 qc -f %s -o allBedpeQc -p %s" % (",".join(allBedpes),
                                                    min(len(allBedpes), cpus))
    callSys([cmd], logger)

    #step 5, unqiue PETs cLoops2 qc
    logger.info("Step5: Unique non-background PETs QC. ")
    cmd = "cLoops2 qc -f %s -o uniNonBgBedpeQc -p %s" % (
        ",".join(uniBedpes), min(len(uniBedpes), cpus))
    callSys([cmd], logger)

    #step 6, combine report
    logger.info("Step5: Generate report. ")
    mata = parseBowtielog()
    matb = pd.read_csv("allBedpeQc_bedpeQc.txt", index_col=0, sep="\t")
    matb.index = [i.split("_all")[0] for i in matb.index]
    matc = pd.read_csv("uniNonBgBedpeQc_bedpeQc.txt", index_col=0, sep="\t")
    matc.index = [i.split("_unique")[0] for i in matc.index]

    for c in matb.columns:
        mata[c] = matb[c]
    mata.to_csv("tracPre_summary.txt", sep="\t")

    mat = {}
    mat["total raw sequences"] = data["totalRaw"]
    mat["after linker removing sequences"] = data["filterLinkers"]
    mat["mapping ratio"] = mata["MappingRatio(%s)"] / 100

    mat["total mapped PETs (mapq>=%s)" % op.mapq] = matb["TotalPETs"]
    mat["total mapped PETs redundancy"] = matb["Redundancy"]
    mat["total mapped PETs intra-chromosomal ratio"] = matb["cisRatio"]
    mat["total mapped PETs close ratio (distance<=1kb)"] = matb["closeRatio"]
    mat["total mapped PETs middle ratio (1kb<distance<=10kb)"] = matb[
        "middleRatio"]
    mat["total mapped PETs distal ratio (10kb<distance)"] = matb["distalRatio"]

    mat["unique noBg PETs"] = matc["TotalPETs"]
    mat["unique noBg mapped PETs intra-chromosomal ratio"] = matc["cisRatio"]
    mat["unique noBg mapped PETs close ratio (distance<=1kb)"] = matc[
        "closeRatio"]
    mat["unique noBg mapped PETs middle ratio (1kb<distance<=10kb)"] = matc[
        "middleRatio"]
    mat["unique noBg mapped PETs distal ratio (10kb<distance)"] = matc[
        "distalRatio"]
    mat = pd.DataFrame(mat)
    mat["Yield"] = mat["unique noBg PETs"].divide(mat["total raw sequences"])

    columns = [
        "total raw sequences",
        "after linker removing sequences",
        "mapping ratio",
        "total mapped PETs (mapq>=%s)" % op.mapq,
        "total mapped PETs redundancy",
        "total mapped PETs intra-chromosomal ratio",
        "total mapped PETs close ratio (distance<=1kb)",
        "total mapped PETs middle ratio (1kb<distance<=10kb)",
        "total mapped PETs distal ratio (10kb<distance)",
        "unique noBg PETs",
        "Yield",
        "unique noBg mapped PETs intra-chromosomal ratio",
        "unique noBg mapped PETs close ratio (distance<=1kb)",
        "unique noBg mapped PETs middle ratio (1kb<distance<=10kb)",
        "unique noBg mapped PETs distal ratio (10kb<distance)",
    ]
    mat = mat[columns]
    mat.to_csv("tracPre_summary.txt", sep="\t")
    cmd = "rm *_bedpeQc.txt"
    os.system(cmd)


if __name__ == '__main__':
    main()
