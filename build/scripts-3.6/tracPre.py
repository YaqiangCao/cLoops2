#!python
#--coding:utf-8--
"""
tracPre.py
Pre-processing code for Trac-looping data, implemented with cLoops2, from fastq to bedpe files and qc report.
2020-02-27: finished and well tested.
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
        help="The directory for raw .fastq.gz files, for example ../1.fastq/ "
    )
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
    parser.add_argument(
        "-mapq",
        dest="mapq",
        required=False,
        default=10,
        type=int,
        help="MAPQ cutoffs for filtering PETs, default is 10."
    )
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
    doBowtie = "bowtie2 -p {cpus} -q --local --very-sensitive -x {ref} -1 {fq1} -2 {fq2} -S {sam}".format(
        cpus=cpus, ref=ref, fq1=fqs[0], fq2=fqs[1], sam=sam)
    logger.info(doBowtie)
    stat, output = subprocess.getstatusoutput(doBowtie)
    #trim with "Warning"
    output = output.split("\n")
    output = [t for t in output if not t.startswith("Warning")]
    output = "\n".join(output)
    logger.info("FLAG_A:" + sample + "\n" + output + "\nFLAG_A\n")
    return sam


def getUniqueBedpe(f, fout):
    """
    Get unique bedpe
    """
    if os.path.isfile(fout):
        return
    print("Getting unique PETs from %s to %s" % (f, fout))
    redus = set()
    #with gzip.open(fout, "wt") as fo:
    with open(fout, "w") as fo:
        #for i, line in enumerate(gzip.open(f, "rt")):
        for i, line in enumerate(open(f)):
            line = line.split("\n")[0].split("\t")
            if len(line) < 6:
                continue
            #remove redudant PETs
            r = hash(tuple(line[:6]))
            if r in redus:
                continue
            else:
                redus.add(r)
            #shroten the name
            #line[6] = str(i)
            fo.write("\t".join(line) + "\n")


def sam2bamBedpe(sam, mapq=10):
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
    return bedpeAll + ".gz"


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

    #mapping
    data = preFqs(op.fqd)
    if len(data) == 0:
        logger.error(
            "No matched _R1.fastq.gz and _R2.fastq.gz in %s. Return." %
            (op.fqd))
        return
    ref = op.ref
    sams = Parallel(n_jobs=op.number,backend="multiprocessing")(
        delayed(tracMapping)(sample, fqs, ref, op.output, cpus=op.cpu)
        for sample, fqs in data.items())
    sams = [sam for sam in sams if sam is not None]

    #sam to bam and bedpe
    cpus = op.number * op.cpu
    ncpus = int(min(len(sams), cpus / 2))
    bedpes = Parallel(n_jobs=ncpus,backend="multiprocessing")(delayed(sam2bamBedpe)(sam) for sam in sams)

    #cLoops2 qc
    cmd = "cLoops2 qc -f %s -o bedpeQc -p %s" % (",".join(bedpes),
                                                 min(len(bedpes), cpus))
    callSys([cmd], logger)

    #combine report
    mata = parseBowtielog()
    matb = pd.read_csv("bedpeQc_bedpeQc.txt", index_col=0, sep="\t")
    matb.index = [i.split("_all")[0] for i in matb.index]
    for c in matb.columns:
        mata[c] = matb[c]
    mata.to_csv("tracPre_summary.txt", sep="\t")
    cmd = "rm bedpeQc_bedpeQc.txt"
    os.system(cmd)


if __name__ == '__main__':
    main()
