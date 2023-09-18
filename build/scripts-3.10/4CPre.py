#!/home/caoy7/anaconda3/envs/astroBoy/bin/python
#--coding:utf-8--
"""
4CPre.py
Pre-processing code for 4C-seq data, implemented with cLoops2, from fastq to fragments and viewpoint bedgraph files.
2022-03-11: finished and well tested.
"""

__author__ = "CAO Yaqiang"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import os
import time
import gzip
import random
import argparse
import subprocess
from glob import glob
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library
import HTSeq
import numpy as np
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
        Preprocess the raw reads of FASTQ files of 4C-seq to reference
        geneome with bowtie2 and obtain the unqiue PETs with quality control
        results.

        Example:
        4CPre.py -fq test -o test -ref ../bowtie2/hg38 -p 5 -mapq 10
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-fq",
                        dest="fq",
                        required=True,
                        type=str,
                        help="The raw .fastq.gz files.")
    parser.add_argument(
        "-o",
        dest="output",
        required=False,
        type=str,
        default="4C",
        help=
        "Output directory, default is 4C, if directory not exists, create one."
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
        "-bait",
        dest="bait",
        required=True,
        type=str,
        help=
        "Bait sequence designed for the view point."
    )
    parser.add_argument(
        "-ligationSite",
        dest="ligationSite",
        required=True,
        type=str,
        help=
        "Ligation site for bait and target sequence. For example, if MboI used, set as GATC."
    )
    parser.add_argument(
        "-genomeFrag",
        dest="genomeFrag",
        required=True,
        type=str,
        help=
        "Genome fragment in bed format. Output of digest_genome.py from HiC-Pro."
    )

    parser.add_argument(
        "-p",
        dest="cpu",
        required=False,
        type=int,
        default=5,
        help="How many cpus used by Bowtie2, default is 5."
    )
    parser.add_argument(
        "-mapq",
        dest="mapq",
        required=False,
        default=10,
        type=int,
        help="MAPQ cutoffs for filtering mapped reads, default is 10."
    )
    parser.add_argument(
        "-cis",
        dest="cis",
        required=False,
        default=False,
        action="store_true",
        help=
        "Whether to only keep intra-chromosomal reads with the bait. The\n"\
        "default is to keep all. "
    )
    parser.add_argument(
        "-log",
        dest="log",
        required=False,
        default=False,
        action="store_true",
        help=
        "Whether to log2 transform the bedGraph signal. Set this to do log2."
    )
    op = parser.parse_args()
    return op


def getBaitPos(bait, ref):
    """
    Get the genomic position for the bait sequence. 
    """
    n = str(random.random())+".fa"
    with open(n,"w") as fo:
        fo.write(">bait\n"+bait+"\n") 
    doBowtie = "bowtie2 --quiet --no-head --no-sq -f --end-to-end -x {ref} -U {fa}".format(ref=ref, fa=n)
    status, output = subprocess.getstatusoutput(doBowtie)
    os.system("rm %s"%n)
    output = output.split("\n")[0].split("\t")
    if output[1] == "16":
        strand = "-"
    else:
        strand = "+"
    chrom = output[2]
    pos = output[3]
    return chrom, pos, strand


def match(sa, sb, miss=2):
    s = 0
    for i in range(len(sa)):
        if sa[i] != sb[i]:
            s += 1
    if s > miss:
        return False
    else:
        return True


def parseSeq(fin, fo, bait, enz, miss=2, rlen=10):
    tot = 0
    q = 0
    with gzip.open(fo, "wt") as fout:
        with gzip.open(fin, "rt") as f:
            for r in FastqGeneralIterator(f):
                r = list(r)
                tot += 1
                #if tot % 100000 == 0:
                #    print("%s reads processed for %s" % (tot, fo))
                s = r[1][:len(bait)]
                m = match(bait, s, miss)
                if m == False:
                    continue
                flag = False
                for i in range(len(bait), len(r[1])):
                    if r[1][i:i + len(enz)] == enz:
                        pos = i + len(enz) + 1
                        if len(r[1]) - pos > rlen:
                            flag = True
                            break
                if flag == False:
                    continue
                r[1] = r[1][pos:]
                r[2] = r[2][pos:]
                q += 1
                fout.write("@%s\n%s\n+\n%s\n" % (r[0], r[1], r[2]))
    return tot, q


def sam2bam(sam, bam):
    """
    SAM to BAM file 
    """
    samview = "samtools view -S %s -b -o %s" % (sam, bam)
    samsort = "samtools sort -@ 2 {bam} -T {pre} -o {bam}".format(
        bam=bam, pre=bam.replace(".bam", ""))
    samindex = "samtools index {bam} {bai}".format(bam=bam,
                                                   bai=bam.replace(
                                                       ".bam", ".bai"))
    rmsam = "rm %s" % (sam)
    cmds = [samview, samsort, samindex, rmsam]
    callSys(cmds, logger)


def doMap(fq, ref,sam,bam, cpus=5):
    #doBowtie = "bowtie2 --no-mixed --no-discordant -p {cpus} -q --local --very-sensitive -x {ref} {fq} -S {sam}".format(
    doBowtie = "bowtie2 -p {cpus} -q --end-to-end --very-sensitive -x {ref} {fq} -S {sam}".format(
        cpus=cpus, ref=ref, fq=fq, sam=sam)
    logger.info(doBowtie)
    status, output = subprocess.getstatusoutput(doBowtie)
    #trim with "Warning"
    output = output.split("\n")
    output = [t for t in output if not t.startswith("Warning")]
    mapRatio = float(output[-1].split("%")[0])
    sam2bam(sam, bam)
    return mapRatio


def bam2Bed(bam, bed, mapq=10):
    """
    Converting BAM file to BED file. 
    bam: bam file path
    bed: bed file path
    mapq: mapq cutoff to remove bad qulity reads.
    """
    fd = os.path.splitext(bed)[0]
    d = os.path.dirname(bed)
    if not os.path.exists(d):
        os.mkdir(d)
    nb = bam.split("/")[-1]
    tmpbam = fd + ".2.bam"
    #important for paired end reads, do it all for all kinds of files.
    samsort = "samtools sort -n -@ 2 {bam} -T {pre} -o {tmpbam}".format(
        bam=bam, tmpbam=nb, pre=nb.replace(".bam", ""))
    rmunmaped = "samtools view -b -q {} -F 4 {} >> {}".format(mapq, nb, tmpbam)
    callSys([samsort, rmunmaped], logger)
    bam2bed = "bamToBed -i {bam} > {bed}".format(bam=tmpbam, bed=bed)
    logger.info(bam2bed)
    status, output = subprocess.getstatusoutput(bam2bed)
    rmbam = "rm {} {}".format(tmpbam, nb)
    callSys([rmbam, "gzip %s" % bed], logger)


def getUniqueBed(f, fout, chrom, cis=False):
    redus = set()
    tot = 0
    c = 0
    with gzip.open(fout, "wt") as fo:
        for line in gzip.open(f, "rt"):
            tot += 1
            line = line.split("\n")[0].split("\t")
            if cis and chrom is not None:
                if line[0] != chrom:
                    continue
            s = int(line[1])
            e = int(line[2])
            r = (line[0], s, e)
            if r in redus:
                continue
            else:
                if line[0] == chrom:
                    c += 1
                redus.add(r)
            fo.write("\t".join(line) + "\n")
    return tot, len(redus), c


def bed2hicFrag(bed, hicFrag, fo, chrom=None, cis=False):
    frags = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for line in open(hicFrag):
        line = line.split("\n")[0].split("\t")
        if cis and chrom is not None and line[0] != chrom:
            continue
        iv = HTSeq.GenomicInterval(line[0], int(line[1]), int(line[2]))
        name = (line[0], line[1], line[2], line[3])
        frags[iv] += name
    c = 0
    with gzip.open(fo, "wt") as fout:
        for line in gzip.open(bed, "rt"):
            line = line.split("\n")[0].split("\t")
            if cis and chrom is not None and line[0] != chrom:
                continue
            strand = line[5]
            if strand == "+":
                p = HTSeq.GenomicPosition(line[0], int(line[1]))
            else:
                p = HTSeq.GenomicPosition(line[0], int(line[2]))
            t = list(frags[p])
            if len(t) == 0:
                continue
            c += 1
            line = t[0]
            fout.write("\t".join(line) + "\n")
    return c


def bed2bdg(f, fout, log=False):
    model = HTSeq.GenomicArray("auto", stranded=False)
    t = 0
    for line in gzip.open(f, "rt"):
        line = line.split("\n")[0].split("\t")
        chrom = line[0]
        s = int(line[1])
        e = int(line[2])
        iv = HTSeq.GenomicInterval(chrom, s, e)
        model[iv] += 1
        t += 1
    with open(fout, "w") as fo:
        for iv, value in model.steps():
            if value > 0:
                value = value / 1.0 / t * 10**6  #RPM
                if log:
                    value = np.log2(value + 1)
                line = [iv.chrom, iv.start, iv.end, value]
                line = list(map(str, line))
                fo.write("\t".join(line) + "\n")


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
    if not os.path.exists(op.fq):
        logger.error("Input %s not exists! Return." % op.fq)
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
            return

    
    logger.info("%s: Start the analysis of sample."%(op.output))
    bait = op.bait.upper()
    enz = op.ligationSite.upper()

    #step 1, get the bait genomic coordinate
    logger.info("%s_Step1: Get bait sequence genomic location."%op.output)
    vpChrom, vpPos, vpStrand = getBaitPos(bait, op.ref)

    #step 2, pre-process fastq files to only keep the reads there are bait and remove the bait sequence
    logger.info("%s_Step2: Trim bait sequence and only keep the target reads."%op.output)
    fastq = op.output+"/"+op.output+".fastq.gz"
    tot,hasBait= parseSeq(op.fq,fastq,bait, enz)
    
    #step 3, mapping the target reads to the genome
    logger.info("%s_Step3: Map the target reads to the reference genome."%op.output)
    sam = op.output+"/"+op.output+".sam"
    bam = op.output+"/"+op.output+".bam"
    mapRatio = doMap(fastq,op.ref,sam,bam,cpus=op.cpu)

    #step 4, get the high quality unqiue reads
    logger.info("%s_Step4: Get the high quality unique reads."%op.output)
    bed = op.output+"/"+op.output+".bed"
    bam2Bed(bam,bed,mapq=op.mapq)
    uniqueBed = op.output+"/"+op.output+"_unique.bed.gz"
    totMapped, uniqueMapped, uniqueCis = getUniqueBed(bed + ".gz", uniqueBed,vpChrom,cis=op.cis)

    #step 5, map the reads to fragments 
    logger.info("%s_Step5: Map the reads to genomic fragments digested."%op.output)
    frag = op.output+"/"+op.output+"_frag.bed.gz" 
    cFrags = bed2hicFrag(uniqueBed,op.genomeFrag,frag,chrom=vpChrom,cis=op.cis)

    #step 6, generate view point bedgraph
    logger.info("%s_Step6: Generate visualization bedGraph file."%op.output)
    bdg = op.output+"/"+op.output+"_frag.bdg"
    bed2bdg(frag, bdg, log=op.log)

    #step 7, generate the qc report
    rs = {
            "0_totalRawReads":tot,
            "1_rawReadsHasBait": hasBait,
            "2_baitRatio": hasBait/tot,
            "3_trimedReadsMappingRatio": mapRatio,
            "4_highQualityMappedReads(MAPQ>=10)":totMapped, 
            "5_highQualityUniqueReads":uniqueMapped,
            "6_redundancy": 1- uniqueMapped/totMapped,
            "7_highQualityUniqueCisReads":uniqueCis,
            "8_cisRatio": uniqueCis/uniqueMapped,
            "9_validFragments": cFrags,
            "10_validRatio": cFrags/uniqueMapped,
        }
    rs = pd.Series(rs)
    rs.to_csv(op.output+"/"+op.output+"_report.txt",sep="\t",header=None)

    logger.info("%s:The analysis finished."%(op.output))


if __name__ == '__main__':
    main()
