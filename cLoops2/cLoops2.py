#!/usr/bin/env python3
#--coding:utf-8 --
"""
cLoops2 main interface to kinds of calling and analysis programme. 
2020-02-18: integrate statBedpe and plotMatrix functions to main part as qc and plot
2020-02-21: integrate aggregrate analysis also in main functions
2020-02-24: integrate loops quantifications in main functions. Change .ixy file to only contain (x,y) information.
2020-02-25: call domains based on correlation matrix. integrate file conversion as dump funciton in main functiuon.
2020-03-05: change call-cis-loops to callLoops, with -trans option for calling intra-chromosomal loops; also change the function name of call-diff-cis-loops to callDiffLoops
2020-03-08: quantification module integrated.
2020-04-06: going to integrate fraction of reads in features module.
2020-05-04: dump to UCSC bigInteract added.
2020-06-08: going to add comp and complot
2020-06-24: going to integrate similarities , comb for combine multiple data sets
2020-08-10: going to add circos plot as Montage module
2020-10-27: integrated interaction density and distance estimation
2020-12-08: integrated find target genes for a set regions/loops through the enhancer promoter network.  
2021-07-16: initial public available for GitHub and PyPI.
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"
__version__ = "0.0.2"

#sys library
import warnings
warnings.filterwarnings("ignore")
import os
import sys
import argparse
from glob import glob
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd
import joblib
import numpy as np
from joblib import Parallel, delayed

#cLoops2, import sorted by module name length
from cLoops2.qc import qcBedpes
from cLoops2.io import parseBedpe, txt2ixy, updateJson, writeNewJson, parseTxt2Loops, combineDirs  #pre-process input BEDPE
from cLoops2.est import estRes, estSat  #estimate reasonable resolution, sequencing depth
from cLoops2.agg import aggPeaks, aggLoops, aggViewPoints, aggDomains, aggTwoAnchors  #aggreate analysis of peaks,loops, viewPoints, domains
from cLoops2.est import getXyDis, getGmmLabelsEps, getKDis, getKDisKneeEps  #estimate eps
from cLoops2.dump import ixy2bed,ixy2bedpe,ixy2hic,ixy2washU,ixy2ucsc,ixy2bdg,ixy2mat #dump files to others
from cLoops2.plot import plotGmmEst, plotKDis, plotKDisE, plotMatHeatmap,plotPETsScatter,plotPETsArches,plotProfiles  #plot
from cLoops2.utils import getLogger  #logger and other utilities
from cLoops2.quant import quantPeaks, quantLoops, quantDomains #quantification of features
from cLoops2.estSim import estSim #estimate similarities 
from cLoops2.estDis import estDis #estimate interaction densities with distance against background 
from cLoops2.filter import filterPETsByPeaks, filterPETsByLoops, filterPETsBySingletons, filterPETsByKNNs, samplePETs  #filter PETs
from cLoops2.montage import montage #montage analysis for rehoboam plot
from cLoops2.callPeaks import callPeaks  #call peaks
from cLoops2.findTargets import findTargets #find regions/loops targets through network
from cLoops2.callDomains import callDomains #call domains based on segregation score
from cLoops2.ano import anaLoops #analysis of loops
from cLoops2.callCisLoops import callCisLoops  #call intra-chromosomal loops
from cLoops2.callDiffLoops import callDiffLoops  #call differential enriched loops between conditions
from cLoops2.callTransLoops import callTransLoops  #call inter-chromosomal loops

#glob settings
#logger
logger = getLogger(os.path.join(os.getcwd(), "cLoops2.log"))
#epilog for argparse
EPILOG = """
Bug reports are welcome and can be put as issue at github repo or sent to 
caoyaqiang0410@gmail.com or yaqiang.cao@nih.gov. Thank you.
"""


def parseEps(eps,emPair=False):
    """
    Parsing command line input eps into order eps.
    @param eps: str
    """
    if "," in str(eps):
        eps = list(map(int, eps.split(",")))
    else:
        eps = [int(eps)]
    if emPair==False:
        eps.sort()  #ascending
    return eps


def parseMinpts(minPts,emPair=False):
    """
    Parsing command line input minPts into order minPts.
    @param minPts: str
    """
    if "," in str(minPts):
        minPts = list(map(int, minPts.split(",")))
    else:
        minPts = [int(minPts)]
    if emPair == False:
        minPts.sort(reverse=True)  #descending
    return minPts


def mainHelp():
    """
    Create the command line interface of the main programme of cLoops2.
    """
    description = """
An enhanced, accurate and flexible peak/domain/loop-calling and analysis tool 
for 3D genomic interaction data.

Use cLoops2 sub-command -h to see detail options and examples for sub-commands.
Available sub-commands are: 
    qc: quality control of BEDPE files before analysis.
    pre: preprocess input BEDPE files into cLoops2 data.
    update: update cLoops2 data files locations.
    combine: combine multiple cLooops2 data directories.
    dump: convert cLoops2 data files to others (BEDPE, HIC, washU, bedGraph and
          contact matrix)
    estEps: estimate eps using Gaussian mixture models or k-distance plot.
    estRes: estimate reasonable contact matrix resolution based on signal 
            enrichment.
    estDis: estimate significant interactions distance range.
    estSat: estimate sequencing saturation based on contact matrix.
    estSim: estimate similarities among samples based on contact matrix.
    filterPETs: filter PETs based on peaks, loops, singleton mode or knn mode. 
    samplePETs: sample PETs according to specific target size.
    callPeaks: call peaks for ChIP-seq, ATAC-seq, ChIC-seq and CUT&Tag or the 
               3D genomic data such as Trac-looping, Hi-TrAC, HiChIP and more.
    callLoops: call loops for 3D genomic data.
    callDiffLoops: call differentially enriched loops for two datasets. 
    callDomains: call domains for 3D genomic data. 
    plot: plot the interaction matrix, genes, view point plot, 1D tracks, 
          peaks, loops and domains for a specific region. 
    montage: analysis of specific regions, producing Westworld Season 3 -like 
             Rehoboam plot. 
    agg: aggregated feature analysis and plots, features can be peaks, view 
         points, loops and domains.
    quant: quantify peaks, loops and domains.
    anaLoops: anotate loops for target genes.
    findTargets: find target genes of genomic regions through networks from 
                 anaLoops.

Examples:
    cLoops2 qc -f trac_rep1.bedpe.gz,trac_rep2.bedpe,trac_rep3.bedpe.gz \\
               -o trac_stat -p 3
    cLoops2 pre -f ../test_GM12878_chr21_trac.bedpe -o trac
    cLoops2 update -d ./trac
    cLoops2 combine -ds ./trac1,./trac2,./trac3 -o trac_combined -keep 1
    cLoops2 dump -d ./trac -o trac -hic
    cLoops2 estEps -d trac -o trac_estEps_gmm -p 10 -method gmm
    cLoops2 estRes -d trac -o trac_estRes -p 10 -bs 25000,5000,1000,200
    cLoops2 estDis -d trac -o trac -plot -bs 1000 
    cLoops2 estSim -ds Trac1,Trac2 -o trac_sim -p 10 -bs 2000 -m pcc -plot
    cLoops2 filterPETs -d trac -peaks trac_peaks.bed -o trac_peaksFiltered -p 10
    cLoops2 samplePETs -d trac -o trac_sampled -t 5000000 -p 10
    cLoops2 callPeaks -d H3K4me3_ChIC -bgd IgG_ChIC -o H3K4me3_cLoops2 -eps 150 \\
                      -minPts 10
    cLoops2 callLoops -d Trac -eps 200,500,1000 -minPts 3 -filter -o Trac -w -j \\
                      -cut 2000
    cLoops2 callLoops -d HiC -eps 1000,5000,10000 -minPts 10,20,50,100 -w -j \\
                      -trans -o HiC_trans 
    cLoops2 callDiffLoops -tloop target_loop.txt -cloop control_loop.txt \\
                          -td ./target -cd ./control -o target_diff
    cLoops2 callDomains -d trac -o trac -bs 10000 -ws 200000
    cLoops2 plot -f test/chr21-chr21.ixy -o test -bs 500 -start 34840000 \\
                 -end 34895000 -triu -1D -loop test_loops.txt -log \\
                 -gtf hg38.gtf -bws ctcf.bw -beds enhancer.bed
    cLoops2 montage -f test/chr21-chr21.ixy -o test -bed test.bed
    cLoops2 agg -d trac -loops trac.loop -peaks trac_peaks.bed \\
                -domains hic_domains.bed -bws CTCF.bw,ATAC.bw -p 20 -o trac 
    cLoops2 quant -d trac -peaks trac_peaks.bed -loops trac.loop \\
                  -domains trac_domain.txt -p 20 -o trac
    cLoops2 anaLoops -loops test_loop.txt -gtf gene.gtf -net -o test
    cLoops2 findTargets -net test_ep_net.sif -tg test_targets.txt \\
                        -bed GWAS.bed -o test 
    More usages and examples are shown when run with cLoops2 sub-command -h.
    """

    #top level, commonly used options
    parser = argparse.ArgumentParser(
        description=description,
        epilog=EPILOG,
        formatter_class=RawTextHelpFormatter,
        prog=None,
        usage=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-d",
        dest="predir",
        required=False,
        default="",
        type=str,
        help=
        "Assign data directory generated by cLoops2 pre to carry out analysis. "
    )
    parser.add_argument(
        "-o",
        dest="fnOut",
        required=False,
        type=str,
        default="cLoops2_output",
        help=
        "Output data directory / file name prefix, default is cLoops2_output.",
    )
    parser.add_argument(
        "-p",
        dest="cpu",
        required=False,
        default=1,
        type=int,
        help=
        "CPUs used to run the job, default is 1, set -1 to use all CPUs\n"\
        "available. Too many CPU could cause out-of-memory problem if there are\n"\
        "too many PETs."
    )
    parser.add_argument(
        "-cut",
        dest="cut",
        required=False,
        default=0,
        type=int,
        help=
        "Distance cutoff to filter cis PETs, only keep PETs with distance\n"\
        ">=cut. Default is 0, no filtering."
    )
    parser.add_argument(
        "-mcut",
        dest="mcut",
        required=False,
        type=int,
        default=-1,
        help= "Keep the PETs with distance <=mcut. Default is -1, no filtering."
    )
    parser.add_argument(
        "-v",
        dest="version",
        action="version",
        version="cLoops2 v%s" % __version__,
        help = "Show cLoops2 verison number and exit."
    )
    #just used to seperate general options and specific options
    parser.add_argument(
        "---",
        dest="version",
        action="version",
        help=
        "Following are sub-commands specific options. This option just show\n"\
        "version of cLoops2.",
        version="cLoops v%s" % __version__,
    )

    #main sub parser
    subparsers = parser.add_subparsers(help=argparse.SUPPRESS)

    #qc
    #quality control of BEDPE files.
    qcDes = """
Get the basic quality control statistical information from interaction BEDPE
files.

Example: 
    cLoops2 qc -f trac_rep1.bedpe.gz,trac_rep2.bedpe,trac_rep3.bedpe.gz -p 3 \\
               -o trac_stat
    """
    qc = subparsers.add_parser(
        'qc',
        description=qcDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )
    qc.add_argument(
        "-f",
        dest="fnIn",
        required=True,
        type=str,
        help=
        "Input BEDPE file(s), .bedpe and .bedpe.gz are both suitable. Multiple\n"\
        "samples can be assigned as -f A.bedpe.gz,B.bedpe.gz,C.bedpe.gz."
    )

    #pre
    #pre-processing of BEDPE files into easy accessed data.
    preDes = """
Preprocess BEDPE PETs into cLoops2 data files.

The output directory contains one .json file for the basic statistics of PETs 
information and .ixy files which are coordinates for every PET. The coordinate
files will be used to call peaks, loops or any other analyses implemented in 
cLoops2. For data backup/sharing purposes, the directory can be saved as 
.tar.gz file through tar. If changed and moved location, run 
***cLoops2 update -d*** to update.

Examples:
    1. keep high quality PETs of chromosome chr21
        cLoops2 pre -f trac_rep1.bepee.gz,trac_rep2.bedpe.gz -o trac -c chr21

    2. keep all cis PETs that have distance > 1kb
        cLoops2 pre -f trac_rep1.bedpe.gz,trac_rep2.bedpe.gz -o trac -mapq 0
    """
    pre = subparsers.add_parser(
        'pre',
        description=preDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )
    pre.add_argument(
        "-f",
        dest="fnIn",
        required=True,
        type=str,
        help=
        "Input BEDPE file(s), .bedpe and .bedpe.gz are both suitable.\n"\
        "Replicates or multiple samples can be assigned as -f A.bedpe.gz,\n"\
        "B.bedpe.gz,C.bedpe.gz to get merged PETs."
    )
    pre.add_argument(
        "-c",
        dest="chroms",
        required=False,
        default="",
        type=str,
        help=
        "Argument to process limited set of chromosomes, specify it as chr1,\n"\
        "chr2,chr3. Use this option to filter reads from such as\n"\
        "chr22_KI270876v1. The default setting is to use the entire set of\n"\
        "chromosomes from the data."
    )
    pre.add_argument(
        "-trans",
        dest="trans",
        required=False,
        default=False,
        action="store_true",
        help=
        "Whether to parse trans- (inter-chromosomal) PETs. The default is to\n"
        "ignore trans-PETs. Set this flag to pre-process all PETs."
    )
    pre.add_argument(
        "-mapq",
        dest="mapq",
        required=False,
        type=int,
        default=10,
        help="MAPQ cutoff to filter raw PETs, default is >=10."
    )

    #update
    updateDes = """
Update cLoops2 data files generated by **cLoops2 pre**.

In the **cLoops2 pre** output directory, there is a .json file annotated with 
the .ixy **absolute paths** and other information. So if the directory is 
moved, or some .ixy files are removed or changed, this command is needed to 
update the paths, otherwise the other analysis modules will not work.

Example:
    cLoops2 update -d ./Trac
    """
    update = subparsers.add_parser(
        'update',
        description=updateDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )

    #combine
    combineDes = """
Combine mulitple cLoops2 data directories.

Example:
    cLoops2 combine -ds ./trac1,./trac2,./trac3 -o combined_trac -keep 1
    """
    combine = subparsers.add_parser(
        'combine',
        description=combineDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )
    combine.add_argument(
        "-ds",
        dest="predirs",
        required=False,
        type=str,
        help=
        "Data directories generated by cLoops2 pre to be combined, seperated by\n"\
        "comma."
    )
    combine.add_argument(
        "-keep",
        dest="keep",
        required=False,
        type=int,
        default=1,
        help=
        "It controls how many PETs will be kept for the same coordinates. If an\n"\
        "integer is given, at most this number of PETs will be kept at the same\n"\
        "location. The default is to keep one PET at the same location. Set 0 to\n"\
        "keep all."
    )

    #dump
    dumpDes = """
Convert cLoops2 data files to other types. Currently supports BED file,BEDPE 
file, HIC file, washU long-range track, bedGraph file and matrix txt file. 

Converting cLoops2 data to .hic file needs "juicer_tools pre" in the command
line enviroment. 
Converting cLoops2 data to legacy washU browser long-range track needs bgzip
and tabix. Format reference: http://wiki.wubrowse.org/Long-range. 
Converting cLoops2 data to UCSC bigInteract track needs bedToBigBed. Format 
reference: https://genome.ucsc.edu/goldenPath/help/interact.html.
Converting cLoops2 data to bedGraph track will normalize value as RPM 
(reads per million). Run with -bdg_pe flag for 1D data such as ChIC-seq,
ChIP-seq and ATAC-seq. 
Converting cLoops2 data to matrix txt file will need specific resolution. 
The output txt file can be loaded in TreeView for visualization or further
analysis. 

Examples:
    1. convert cLoops2 data to single-end .bed file fo usage of BEDtools or 
       MACS2 for peak-calling with close PETs
        cLoops2 dump -d trac -o trac -bed -mcut 1000

    2. convert cLoops2 data to .bedpe file for usage of BEDtools, only keep 
       PETs distance >1kb and < 1Mb
        cLoops2 dump -d trac -o trac -bedpe -bedpe_ext -cut 1000 -mcut 1000000 

    3. convert cLoops2 data to .hic file to load in juicebox
        cLoops2 dump -d trac -o trac -hic -hic_org hg38 \\
                    -hic_res 200000,20000,5000
    
    4. convert cLoops2 data to washU long-range track file, only keep PETs 
       distance > 1kb 
        cLoops2 dump -d trac -o trac -washU -washU_ext 50 -cut 1000
    
    5. convert cLoops2 data to UCSC bigInteract track file 
        cLoops2 dump -d trac -o trac -ucsc -ucsc_cs ./hg38.chrom.sizes 

    6. convert interacting cLoops2 data to bedGraph file with all PETs
        cLoops2 dump -d trac -o trac -bdg -bdg_ext 100

    7. convert 1D cLoops2 data (such as ChIC-seq/ChIP-seq/ATAC-seq) to bedGraph 
       file 
        cLoops2 dump -d trac -o trac -bdg -pe 

    8. convert 3D cLoops2 data (such as Trac-looping) to bedGraph file for peaks
        cLoops2 dump -d trac -o trac -bdg -mcut 1000

    9. convert one region in chr21 to contact matrix correlation matrix txt file 
        cLoops2 dump -d test -mat -o test -mat_res 10000 \\
                    -mat_chrom chr21-chr21 -mat_start 36000000 \\
                    -mat_end 40000000 -log -corr
    """
    dump = subparsers.add_parser(
        'dump',
        description=dumpDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )
    dump.add_argument(
        "-bed",
        dest="bed",
        required=False,
        action="store_true",
        default=False,
        help= "Convert data to single-end BED file."
    )
    dump.add_argument(
        "-bed_ext",
        dest="bed_ext",
        required=False,
        type=int,
        default=50,
        help=
        "Extension from the center of the read to both ends for BED file.\n"\
        "Default is 50."
    )
    dump.add_argument(
        "-bedpe",
        dest="bedpe",
        required=False,
        action="store_true",
        default=False,
        help= "Convert data to BEDPE file."
    )
    dump.add_argument(
        "-bedpe_ext",
        dest="bedpe_ext",
        required=False,
        type=int,
        default=50,
        help=
        "Extension from the center of the PET to both ends for BEDPE file.\n"\
        "Default is 50."
    )
    dump.add_argument(
        "-hic",
        dest="hic",
        required=False,
        action="store_true",
        default=False,
        help= "Convert data to .hic file."
    )
    dump.add_argument(
        "-hic_org",
        dest="hic_org",
        required=False,
        type=str,
        default="hg38",
        help=
        "Organism required to generate .hic file,default is hg38. If the\n"\
        "organism is not available, assign a chrom.size file."
    )
    dump.add_argument(
        "-hic_res",
        dest="hic_res",
        type=str,
        default="1000,5000,25000,50000,100000,200000",
        help= 
        "Resolutions used to generate .hic file. Default is 1000,5000,25000,\n"\
        "50000,100000,200000."
    )
    dump.add_argument(
        "-washU",
        dest="washU",
        required=False,
        action="store_true",
        default=False,
        help= 
        "Convert data to legacy washU browser long-range track."
    )
    dump.add_argument(
        "-washU_ext",
        dest="washU_ext",
        required=False,
        type=int,
        default=50,
        help=
        "Extension from the center of the PET to both ends for washU track.\n"\
        "Default is 50."
    )
    dump.add_argument(
        "-ucsc",
        dest="ucsc",
        required=False,
        action="store_true",
        default=False,
        help= "Convert data to UCSC bigInteract file track."
    )
    dump.add_argument(
        "-ucsc_ext",
        dest="ucsc_ext",
        required=False,
        type=int,
        default=50,
        help="Extension from the center of the PET to both ends for ucsc\n"\
        "track. Default is 50."
    )
    dump.add_argument(
        "-ucsc_cs",
        dest="ucsc_cs",
        required=False,
        default="",
        type=str,
        help="A chrom.sizes file. Can be obtained through fetchChromSizese.\n"\
        "Required for -ucsc option."
    )
    dump.add_argument(
        "-bdg",
        dest="bdg",
        required=False,
        action="store_true",
        default=False,
        help= "Convert data to 1D bedGraph track file."
    )
    dump.add_argument(
        "-bdg_ext",
        dest="bdg_ext",
        required=False,
        type=int,
        default=50,
        help="Extension from the center of the PET to both ends for\n"\
        "bedGraph track. Default is 50."
    )
    dump.add_argument(
        "-bdg_pe",
        dest="bdg_pe",
        required=False,
        action="store_true",
        default=False,
        help= 
        "When converting to bedGraph, argument determines whether to treat PETs\n"\
        "as ChIP-seq, ChIC-seq or ATAC-seq paired-end libraries. Default is not.\n"
        "PETs are treated as single-end library for interacting data."
    )
    dump.add_argument(
        "-mat",
        dest="mat",
        required=False,
        action="store_true",
        default=False,
        help= 
        "Convert data to matrix txt file with required resolution."
    )
    dump.add_argument(
        "-mat_res",
        dest="mat_res",
        required=False,
        default=5000,
        type=int,
        help=
        "Bin size/matrix resolution (bp) to generate the contact matrix. \n"\
        "Default is 5000 bp. "
    )
    dump.add_argument(
        "-mat_chrom",
        dest="chrom",
        required=False,
        default="",
        type=str,
        help=
        "The chrom-chrom set will be processed. Specify it as chr1-chr1.\n"\
    )
    dump.add_argument(
        "-mat_start",
        dest="start",
        required=False,
        type=int,
        default=-1,
        help=
        "Start genomic coordinate for the target region. Default will be the\n"\
        "smallest coordinate from specified chrom-chrom set."
    )
    dump.add_argument(
        "-mat_end",
        dest="end",
        required=False,
        type=int,
        default=-1,
        help=
        "End genomic coordinate for the target region. Default will be the\n"\
        "largest coordinate from specified chrom-chrom set."
    )
    dump.add_argument(
        "-log",
        dest="log",
        required=False,
        action="store_true",
        default=False,
        help=
        "Whether to log transform the matrix. Default is not."
    )
    dump.add_argument(
        "-m",
        dest="method",
        type=str,
        choices=["obs", "obs/exp"],
        default="obs",
        help=
        "The type of matrix, observed matrix or observed/expected matrix, \n"\
        "expected matrix will be generated by shuffling PETs. Default is\n"\
        "observed."
    )
    dump.add_argument(
        "-corr",
        dest="corr",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to get the correlation matrix. Default is not. "
    )
    dump.add_argument(
        "-norm",
        dest="norm",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to normalize the matrix with z-score. Default is not."
    )

    #estEps
    #estimation eps using estimated peak size GMM or k-distance plot
    estEpsDes = """
Estimate key parameter eps. 

Two methods are implemented: 1) unsupervised Gaussian mixture model (gmm), and 
2) k-distance plot (k-dis,-k needed). Gmm is based on the assumption that PETs 
can be classified into self-ligation (peaks) and inter-ligation (loops). K-dis
is based on the k-nearest neighbors distance distribution to find the "knee", 
which is where the distance (eps) between neighbors has a sharp increase along
the k-distance curve. K-dis is the traditional approach literatures, but it is
much more time consuming than gmm, and maybe only fit to small cases. If both 
methods do not give nice plots, please turn to the empirical parameters you 
like, such as 100,200 for ChIP-seq -like data, 5000,1000 for Hi-C and etc.

Examples: 
    1. estimate eps with Gaussian mixture model    
        cLoops2 estEps -d trac -o trac_estEps_gmm -p 10 -method gmm

    2. estimate eps with k-nearest neighbors distance distribution
        cLoops2 estEps -d trac -o trac_estEps_kdis -p 10 -method k-dis -k 5
    """

    estEps = subparsers.add_parser(
        'estEps',
        description=estEpsDes,
        add_help=False,
        formatter_class=RawTextHelpFormatter,
        parents=[parser],
        usage=argparse.SUPPRESS,
    )
    estEps.add_argument(
        "-fixy",
        dest="fixy",
        required=False,
        default="",
        type=str,
        help=
        "Assign the .ixy file to estimate eps inside of the whole directory\n"\
        "generated by cLoops2 pre. For very large data, especially Hi-C, this\n"\
        "option is recommended for chr1 (or the smaller one) to save time."
    )
    estEps.add_argument(
        "-k",
        dest="knn",
        required=False,
        #default=5,
        default=0,
        type=int,
        help=
        "The k-nearest neighbors used to draw the k-distance plot. Default is 0\n"\
        "(not running), set this when -method k-dis. Suggested 5 for\n"\
        "ChIA-PET/Trac-looping data, 20 or 30 for Hi-C like data."
    )
    estEps.add_argument(
        "-method",
        dest="epsMethod",
        required=False,
        default="gmm",
        choices=["gmm", "k-dis"],
        help=
        "Two methods can be chosen to estimate eps. Default is Gmm. See above\n"\
        "for difference of the methods."
    )

    #estimation of reasonable contact matrix resolution based on signal enrichment
    estResDes = """
Estimate reasonable genome-wide contact matrix resolution based on signal 
enrichment. 

PETs will be assigned to contact matrix bins according to input resolution. A 
bin is marked as [nx,ny], and a PET is assigned to a bin by nx = int((x-s)/bs),
ny = int((y-s)/bs), where s is the minimal coordinate for all PETs and bs is 
the bin size. Self-interaction bins (nx=ny) will be ignored. The bins only 
containing singleton PETs are assumed as noise. 

The output is a PDF plot, for each resolution, a line is separated into two 
parts: 1) dash line indicated linear increased trend of singleton PETs/bins; 2)
solid thicker line indicated non-linear increased trend of higher potential 
signal PETs/bins. The higher the ratio of signal PETs/bins, the easier it it to
find loops in that resolution. The closer to the random line, the higher the 
possibility to observe evenly distributed signals.  

We expect the highest resolution with >=50% PETs are not singletons.

Example:
    cLoops2 estRes -d trac -o trac -bs 10000,5000,1000 -p 20
"""
    estRes = subparsers.add_parser(
        'estRes',
        description=estResDes,
        formatter_class=RawTextHelpFormatter,
        add_help=False,
        parents=[parser],
        usage=argparse.SUPPRESS,
    )
    estRes.add_argument(
        "-bs",
        dest="binSize",
        required=True,
        type=str,
        default="5000",
        help=
        "Candidate contact matrix resolution (bin size) to estimate signal\n"\
        "enrichment. A series of comma-separated values or a single value can\n"\
        "be used as input. For example,-bs 1000,5000,10000. Default is 5000."
    )

    #estimation of interaction distance limitaion
    estDisDes = """
Estimate the significant interaction distance limitation by getting the observed
and expected random background of the genomic distance vs interaction frequency.

Example:
    cLoops2 estDis -d trac -o trac -bs 5000 -p 20 -plot
    """
    estDis = subparsers.add_parser(
        'estDis',
        description=estDisDes,
        formatter_class=RawTextHelpFormatter,
        add_help=False,
        parents=[parser],
        usage=argparse.SUPPRESS,
    )
    estDis.add_argument(
        "-c",
        dest="chroms",
        required=False,
        default="",
        type=str,
        help=
        "Whether to process limited chroms, specify it as chr1,chr2,chr3, \n"\
        "default is not. Use this to save time for quite big data."
    )
    estDis.add_argument(
        "-bs",
        dest="binSize",
        required=False,
        default=5000,
        type=int,
        help=
        "Bin size / contact matrix resolution (bp) to generate the contact\n"\
        "matrix for estimation, default is 5000 bp."
    )
    estDis.add_argument(
        '-r',
        dest="repeats",
        required=False,
        default=10,
        type=int,
        help=
        "The reapet times to shuffle PETs to get the mean expected background,\n"\
        "default is 10."
    )
    estDis.add_argument('-plot',
                        dest="plot",
                        required=False,
                        action="store_true",
                        help="Set to plot the result.")

    #estimation of sequencing saturation based on signal enrichment of contact matrix
    estSatDes = """
    Estimate reasonable sequencing saturation based on contact matrix. 

    Example:
        cLoops2 estSat -d trac -o trac -bs 5000 -p 20
    """
    estSat = subparsers.add_parser(
        'estSat',
        description=estSatDes,
        formatter_class=RawTextHelpFormatter,
        add_help=False,
        parents=[parser],
        usage=argparse.SUPPRESS,
    )
    estSat.add_argument(
        "-bs",
        dest="binSize",
        required=True,
        type=str,
        default="5000",
        help=
        "Candidate contact matrix resolution (bin size) to estimate sequencing\n"\
        "saturation. A series of comma-separated values or a single value can\n"\
        "be used as input. For example, -bs 1000,5000,10000. Default is 5000."
    )
    estSat.add_argument(
        "-tol",
        dest="tol",
        required=False,
        type=int,
        default=5,
        help=
        "How many PETs are needed to define a detected contact matrix bin.\n"
        "Default is 5"
    )


    #estimation of similarities 
    estSimDes = """
Estimate the similarity of >=two experiments based on contact matrix.

Two methods are implemented: 1) PCC, the contact matrix are flatten as vectors,
then Pearson correlation coeificient is used to measure the similarities. 2) 
PCA, the contact matrix are first combined and averaged as C, then matrix C is
used for PCA embeding features, each sample contact matrix are projected to the
reduced dimensions of C, which are furthur flatten and used to caculate Pearson
correlation coeificient. PCA based method is much more time and memory 
consuming than direct PCC. 

Example:
    1. run with PCC, two samples will generate a scatter plot
        cLoops2 estSim -ds Trac1,Trac2 -o trac_sim -cut 0 -p 10 -bs 1000 -m pcc \\
                       -plot
    
    2. run with PCA, more than two samples will generate a correlation heatmap
        cLoops2 estSim -ds Trac1,Trac2,Trac3 -o trac_sim -cut 0 -p 10 -bs 2000 \\
                       -m pcc -plot

    3. control noise, run with all PETs 
        cLoops2 estSim -ds Trac1,Trac2 -o trac_sim -p 10 -bs 1000 -pcut 0 -plot
    
    4. control noise, run with bins >= signal PETs, for example 2
        cLoops2 estSim -ds Trac1,Trac2 -o trac_sim -p 10 -bs 1000 -pcut 2 -plot

    5. get the data, and plot by yourself
        cLoops2 estSim -ds Trac1,Trac2 -o trac_sim -p 10 -bs 1000 -pcut 2
"""
    estSim = subparsers.add_parser(
        'estSim',
        description=estSimDes,
        formatter_class=RawTextHelpFormatter,
        add_help=False,
        parents=[parser],
        usage=argparse.SUPPRESS,
    )
    estSim.add_argument(
        "-ds",
        dest="predirs",
        required=False,
        type=str,
        help=
        "Data directories generated by cLoops2 pre to be combined, seperated by\n"\
        "comma."
    )
    estSim.add_argument(
        "-bs",
        dest="binSize",
        required=False,
        type=int,
        default=5000,
        help=
        "Bin size (bp) to generate the contact matrix for estimation. Only one\n"\
        "int value is supported. Default is 5000 bp."
    )
    estSim.add_argument(
        "-m",
        dest="method",
        type=str,
        choices=["pcc", "pca"],
        default="pcc",
        help=
        "Method to compare the data, default is pcc (Pearson Correlation\n"\
        "Coefficient). Another option is pca."
    )
    estSim.add_argument(
        "-pcut",
        dest="pcut",
        type=int,
        default=2,
        required=False,
        help=
        "If the method pcc is used, remove the bins that no sample has PETs >\n"\
        "pcut. Default is 2."
    )
    estSim.add_argument(
        "-n",
        dest="n_components",
        type=int,
        default=2,
        required=False,
        help=
        "The top n_components of PCA embeding features if pca is selected.\n"\
        "Default is 2."
    )
    estSim.add_argument(
        "-plot",
        dest="plot",
        required=False,
        default=False,
        action="store_true",
        help="Whether to plot the similarities. Default is not."
    )

    #filter PETs to reduce noise.
    filterPETsDes = """
Filter PETs according to peaks/domains/loops/singletons/KNNs. 

If any end of the PETs overlap with features such as peaks or loops, the PET 
will be kept. Filtering can be done before or after peak/loop-calling. Input 
can be peaks or loops, but should not be be mixed. The -singleton mode is based
on a specified contact matrix resolution, if there is only one PET in the bin, 
the singleton PETs will be filtered. The -knn is based on noise removing step 
of blockDBSCAN. 

Examples:
    1. keep PETs overlapping with peaks
        cLoops2 filterPETs -d trac -peaks peaks.bed -o trac_filtered

    2. keep PETs that do not overlap with any blacklist regions
        cLoops2 filterPETs -d trac -peaks bg.bed -o trac_filtered -iv

    3. keep PETs that overlap with loop anchors
        cLoops2 filterPETs -d trac -loops test_loops.txt -o trac_filtered

    4. keep PETs that both ends overlap with loop anchors
        cLoops2 filterPETs -d trac -loops test_loops.txt -o trac_filtered -both

    5. keep non-singleton PETs based on 1kb contact matrix
        cLoops2 filterPETs -d trac -o trac_filtered -singleton -bs 1000

    6. filter PETs based on blockDBSCAN knn noise removing
        cLoops2 filterPETs -d trac -o trac_filtered -knn -eps 1000 -minPts 5

"""
    filterPETs = subparsers.add_parser(
        'filterPETs',
        description=filterPETsDes,
        formatter_class=RawTextHelpFormatter,
        parents=[parser],
        add_help=False,
        usage=argparse.SUPPRESS,
    )
    filterPETs.add_argument(
        "-peaks",
        dest="fbed",
        default="",
        required=False,
        help=
        "BED file of genomic features (such as promoters, enhancers, ChIP-seq,\n"\
        "ATAC-seq peaks,TADs) to filter PETs."
    )
    filterPETs.add_argument(
        "-loops",
        dest="floop",
        required=False,
        default="",
        type=str,
        help=
        "The loop.txt file generated by cLoops2, can be loops or domains, to\n"\
        "filter PETs."
    )
    filterPETs.add_argument(
        "-gap",
        dest="gap",
        required=False,
        default=1,
        type=int,
        help=
        "If the distance between two genomic features is <=gap, the two regions\n"
        "will be combined. Default is 1. Set to >=1."
    )
    filterPETs.add_argument(
        "-singleton",
        dest="singleton",
        required=False,
        action="store_true",
        default=False,
        help=
        "Whether to use singleton mode to filter PETs. Contact matrix\n"\
        "resolution with -bs is required. Singleton PETs in contact matrix bins\n"
        "will be filtered."
    )
    filterPETs.add_argument(
        "-bs",
        dest="binSize",
        required=False,
        type=int,
        default=5000,
        help=
        "The contact matrix bin size for -singleton mode filtering. Default is\n"\
        "5000."
    )
    filterPETs.add_argument(
        "-knn",
        dest="knn",
        required=False,
        action="store_true",
        default=False,
        help=
        "Whether to use noise removing method in blockDBSCAN to filter PETs,\n"\
        "-eps and -minPts are required."
    )
    filterPETs.add_argument(
        "-eps",
        dest="eps",
        default=1000,
        required=False,
        type=int,
        help=
        "Same to callPeaks and callLoops, only used to filter PETs for -knn\n"\
        "mode. Default is 1000. Only one value is supported."
    )
    filterPETs.add_argument(
        "-minPts",
        dest="minPts",
        default=5,
        type=int,
        help=
        "Same to callPeaks and callLoops, only used to filter PETs for -knn\n"\
        "mode. Default is 5. Only one value is supported."
    )
    filterPETs.add_argument(
        "-iv",
        dest="iv",
        required=False,
        action="store_true",
        default=False,
        help=
        "Whether to only keep PETs not in the assigned regions, behaves like\n"
        "grep -v."
    )
    filterPETs.add_argument(
        "-both",
        dest="both",
        required=False,
        action="store_true",
        default=False,
        help=
        "Whether to only keep PETs that both ends overlap with loop anchors.\n"\
        "Default is not."
    )

 
    #sample PETs to similar depth for fair comparasion
    samplePETsDes = """
Sampling PETs to target total size. 

If there are multiple sample libraries and the total sequencing depths vary a 
lot, and you want to compare the data fairly, it's better to sample them to 
similar total PETs (either down-sampling or up-sampling), then call peaks/loops
with the same parameters. 

Example:
    cLoops2 samplePETs -d trac -o trac_sampled -tot 5000000 -p 10
    """
    samplePETs = subparsers.add_parser(
        'samplePETs',
        parents=[parser],
        add_help=False,
        usage=argparse.SUPPRESS,
        description=samplePETsDes,
        formatter_class=RawTextHelpFormatter,
    )
    samplePETs.add_argument(
        "-tot",
        dest="tot",
        default=0,
        required=True,
        type=int,
        help="Target total number of PETs.",
    )

    #calling intra-chromosomal loops
    callPeaksDes = """
Call peaks based on clustering. 

Well tested work for ChIP-seq, ChIC-seq, ATAC-seq, CUT&RUN -like or the 3D
genomic data such as Hi-TrAC/Trac-looping, ChIA-PET and HiChIP.

There are three steps in the algorithm: 1) cluster the PETs to find 
self-ligation clusters, which are candidate peaks; 2) estimate the significance
of candidate peaks with local background; 3) if given control data, further 
compare candidate peaks to control data. If running multiple clusterings with
separated parameters, the clusters will be combined and callPeaks will output 
the most significant one based on overlaps. 

Key parameters are -eps and -minPts, both are key parameters in the clustering
algorithm blockDBSCAN. Eps indicates the distance that define two points (PETs) 
being neighbors, while minPts indicatess the minial number of points required 
for a cluster to form.  For sharp-peak like data (ATAC-seq, TF ChIC-seq), set
-eps small such as 100 or 150. For broad-peak like data, such as H3K27me3 
ChIP-seq and ChIC-seq, set -eps large as 500 or 1000. 

Eps will affect more than minPts for sensitivity.

Examples:
    1. call peaks for Trac-looping  
        cLoops2 callPeaks -d trac -eps 100 -minPts 10 -o trac -p 10

    2. call peaks for sharp-peak like ChIC-seq without control data
        cLoops2 callPeaks -d ctcf_chic -o ctcf_chic -p 10

    3. call peaks for broad-peak like ChIC-seq with IgG as control
        cLoops2 callPeaks -d H3K27me3 -bgd IgG -eps 500,1000 -minPts 10 \\
                          -o H3K27me3 

    4. call peaks for sharp-peak ChIC-seq with linear fitting scaled control 
       data
        cLoops2 callPeaks -d ctcf -bgd IgG -eps 150 -minPts 10 -o ctcf -p 10\\
                          -bgm lf

    5. call peaks with sentitive mode to get comprehensive peaks for CUT&TAG
        cLoops2 callPeaks -d H3K27ac -bgd IgG -sen -p 10

    6. filter PETs first and then call peaks for H3K27ac HiChIP, resulting much
       better accurate peaks
        cLoops2 filterPETs -d h3k27ac_hichip -o h3k27ac_hichip_filtered -knn \\
                           -eps 500 -minPts 5
        cLoops2 callPeaks -d h3k27ac_hichip_filtered -eps 200,500 -minPts 10 \\
                          -p 10

    7. call peaks for interaction data as single-end data 
        cLoops2 callPeaks -d h3k27ac -o h3k27ac -split -eps 200,500 -minPts 10 \\
                          -p 10

    8. call differential peaks between WT and KO condition
        cLoops2 callPeaks -d MLL4_WT -bgd MLL4_KO -o MLL4_WTvsKO -p 10
        cLoops2 callPeaks -d MLL4_KO -bgd MLL4_WT -o MLL4_KOvsWT -p 10
    """
    callPeaks = subparsers.add_parser(
        'callPeaks',
        parents=[parser],
        add_help=False,
        usage=argparse.SUPPRESS,
        description=callPeaksDes,
        formatter_class=RawTextHelpFormatter,
    )
    callPeaks.add_argument(
        "-eps",
        dest="eps",
        default="100,200",
        required=False,
        type=str,
        help=
        "Distance that defines two points (PETs) being neighbors, eps in\n"\
        "blockDBSCAN as key parameter, multiple eps can be assigned such as\n"\
        "100,200,300 to run multiple clusterings, the results will be combined.\n"\
        "For callPeaks, the default is 100,200. If the data show much more broad\n"\
        "feature such as H3K27me3 and H3K4me1, increase it to 500,1000 or larger.\n"\
        "If expecting both narrow and broad peaks in the data, set -eps 100,200,\n"\
        "500,1000."
    )
    callPeaks.add_argument(
        "-minPts",
        dest="minPts",
        default="5",
        type=str,
        help=
        "Points required in a cluster, minPts in blockDBSCAN, key parameter,\n"\
        "multiple minPts can be assigned such as 3,5 to run multiple\n"\
        "clusterings, the results will be combined. For callPeaks, the default\n"\
        "is 5. If the data have many reads, increasing minPts such as 10,20."
    )
    callPeaks.add_argument(
        "-pcut",
        dest="pcut",
        default=1e-2,
        type=float,
        help=
        "Bonferroni corrected Poisson p-value cutoff to determine significant\n"\
        "peaks. Default is 1e-2."
    )
    callPeaks.add_argument(
        "-bgd",
        dest="bgd",
        default="",
        type=str,
        help=
        "Assign control data (IgG, Input) directory generated by cLoops2 pre to\n"\
        "carry out analysis. Default is no background.",
    )
    callPeaks.add_argument(
        "-bgm",
        dest="bgm",
        choices=["ratio", "lf"],
        default="lf",
        type=str,
        help=
        "How to scale the target data with control data. Available options are\n"\
        "'ratio' and 'lf'. 'ratio' is based on library size and 'lf' means\n"\
        "linear fitting for control and target candidate peaks nearby regions.\n"\
        "Default is 'lf'. The scaling factor estimated by lf usually is a little\n"\
        "larger than ratio. In other words, the higher the scaling factor, the\n"\
        "less sensitive the results. Default is lf."
    )
    callPeaks.add_argument(
        "-pseudo",
        dest="pseudo",
        default=1,
        type=int,
        help=
        "Pseudo counts for local background or control data to estimate the\n"\
        "significance of peaks if no PETs/reads in the background. Default is\n"\
        "1. Set it larger for noisy data, 0 is recommend for very clean data\n"\
        "such as well prepared CUT&Tag."
    )
    callPeaks.add_argument(
        "-sen",
        dest="sen",
        required=False,
        action="store_true",
        help=
        "Whether to use sensitive mode to call peaks. Default is not. If only a\n"\
        "few peaks were called, while a lot more can be observed\n"\
        "from visualization, try this option. Adjust -pcut or filter by\n"
        "yourself to select significant ones."
    )
    callPeaks.add_argument(
        "-split",
        dest="split",
        required=False,
        action="store_true",
        help=
        "Whether to split paired-end as single end data to call peaks. Sometimes\n"\
        "works well for Hi-TrAC/Trac-looping and HiChIP. Can get more peaks."
    )
    callPeaks.add_argument(
        "-splitExt",
        dest="splitExt",
        required=False,
        type=int,
        default=50,
        help=
        "When run with -split, the extension to upstraem and downstream, \n"\
        "default is 50."
    )

    #calling intra-chromosomal loops
    callLoopsDes = """
Call loops based on clustering. 

Well tested work for Hi-TrAC/TrAC-looping, HiCHiP, ChIA-PET and Hi-C.

Similar to call peaks, there are three main steps in the algorithm: 1) cluster 
the PETs to find inter-ligation clusters, which are candidate loops; 2) 
estimate the significance of candidate loops with permutated local background. 
3) If -hic option not selected, the loop anchors will be checked for peak-like 
features, only peak-like anchors are kept. If running multiple clusterings, 
the clusters will be combined and callLoops will output the most significant 
one based on overlaps. 

Similar to callPeaks, keys parameters are -eps and -minPts. For sharp-peak like 
interaction data, set -eps small such as 500,1000. For broad-peak like data, 
such as H3K27ac HiChIP, set -eps big as 1000,2000. For Hi-C and HiChIP data, 
bigger -minPts is also needed, such as 20,50. 

Please note that the blockDBSCAN implementation in cLoops2 is much more 
sensitive than cDBSCAN in cLoops, so the same parameters can generate quite 
different results. With -hic option, cDBSCAN will be used. 

Examples:
    1. call loops for Hi-TrAC/Trac-looping
        cLoops2 callLoops -d trac -o trac -eps 200,500,1000,2000 -minPts 5 -w -j

    2. call loops for Hi-TrAC/Trac-looping with filtering short distance PETs 
       and using maximal estimated distance cutoff
        cLoops2 callLoops -d trac -o trac -eps 200,500,1000,2000 -minPts 5 \\
                          -cut 1000 -max_cut -w -j

    3. call loops for Hi-TrAC/Trac-looping and get the PETs with any end 
       overlapping loop anchors
        cLoops2 callLoops -d trac -o trac -eps 200,500,1000,2000 -minPts 5 -w \\
                          -j -filterPETs

    4. call loops for high-resolution Hi-C like data 
        cLoops2 callLoops -d hic -o hic -eps 2000,5000,10000 -minPts 20,50 -w -j
    
    5. call inter-chromosomal loops (for most data, there will be no significant 
       inter-chromosomal loops)
        cLoops2 callLoops -d HiC -eps 5000 -minPts 10,20,50,100,200 -w -j -trans\\                          
                          -o HiC_trans
    """
    callLoops = subparsers.add_parser(
        'callLoops',
        parents=[parser],
        add_help=False,
        usage=argparse.SUPPRESS,
        description=callLoopsDes,
        formatter_class=RawTextHelpFormatter,
    )
    callLoops.add_argument(
        "-eps",
        dest="eps",
        default=0,
        required=False,
        help=
        "Distance that defines two points (PETs) being neighbors, eps in\n"\
        "blockDBSCAN as key parameter, multiple eps can be assigned such as\n"\
        "200,500,1000,2000 to run multiple clusterings, the results will be\n"\
        "combined. No default value, please give the input."
    )
    callLoops.add_argument(
        "-minPts",
        dest="minPts",
        default=0,
        help=
        "Points required in a cluster. minPts in blockDBSCAN is a key parameter.\n"
        "Empirically 5 is good for TFs and histone modification ChIA-PET data\n"
        "and Trac-looping. For data like HiChIP and Hi-C, set it larger, like\n"
        ">=20. The input can be a series, and the final loops will have the\n"
        "PETs>= max(minPts). ")
    callLoops.add_argument(
        "-plot",
        dest="plot",
        required=False,
        action="store_true",
        help=
        "Whether to plot estimated inter-ligation and self-ligation PETs\n"\
        "distance distribution. Default is not to generate a plot."
    )
    callLoops.add_argument(
        "-i",
        dest="ucsc",
        required=False,
        action="store_true",
        help=
        "Whether to convert loops to UCSC Interact track to visualize in UCSC.\n"\
        "Default is not, set this flag to save."
    )
    callLoops.add_argument(
        "-j",
        dest="juicebox",
        required=False,
        action="store_true",
        help=
        "Whether to convert loops to 2D feature annotations to visualize in\n"\
        "Juicebox. Default is not, set this flag to save."
    )
    callLoops.add_argument(
        "-w",
        dest="washU",
        required=False,
        action="store_true",
        help=
        "Whether to save tracks of loops to visualize in legacy and new washU.\n"\
        "Default is not, set this flag to save two files."
    )
    callLoops.add_argument(
        "-max_cut",
        dest="max_cut",
        required=False,
        action="store_true",
        help=
        "When running cLoops with multiple eps or minPts, multiple distance\n"\
        "cutoffs for self-ligation and inter-ligation PETs will be estimated\n"\
        "based on the overlaps of anchors. Default option is the minimal one\n"\
        "will be used to filter PETs for candidate loop significance test.\n"\
        "Set this flag to use maximal one, will speed up for significance test.\n"
    )
    callLoops.add_argument(
        "-hic",
        dest="hic",
        required=False,
        action="store_true",
        help=
        "Whether to use statistical cutoffs for Hi-C to output significant loops.\n"\
        "Default is not, set this option to enable. Additionally, with -hic\n"\
        "option, there is no check for anchors requiring they looking like peaks."
    )
    callLoops.add_argument(
        "-filter",
        dest="filterPETs",
        required=False,
        action="store_true",
        help=
        "Whether to filter raw PETs according to called loops. The filtered\n"\
        "PETs can show clear view of interactions or be used to call loops again."
    )
    callLoops.add_argument(
        "-trans",
        dest="trans",
        required=False,
        default=False,
        action="store_true",
        help=
        "Whether to call trans- (inter-chromosomal) loops. Default is not, set\n"\
        "this flag to call. For most common cases, not recommended, only for\n"\
        "data there are obvious visible trans loops."
    )
    callLoops.add_argument(
        "-emPair",
        dest="emPair",
        required=False,
        action="store_true",
        help=
        "By default eps and minPts combinations will be used to run clustering.\n"\
        "With this option, for example eps=500,1000 and minPts=5,10, only (500,5)\n"\
        "and (1000,10) as parameters of clustering will be run. Input number of\n"\
        "eps and minPts should be same."
    )
 
    callDiffLoopsDes = """
Call differentially enriched intra-chromosomal loops between two conditions.

Similar to calling peaks with control data, calling differentially enriched 
loops is based on scaled PETs and the Poisson test. There are three main steps 
in the algorithm: 1) merge the overlapped loops, quantify them and their 
permutated local background regions; 2) fit the linear transformation of 
background target interaction density to control background data based on 
MANorm2; therefore, if there are more than than two samples, others can be 
scaled to the reference sample for quantitative comparison; 3) estimate the 
fold change (M) cutoff and average (A) cutoff using the background data with 
the control of FDR, assuming there should no differentially significant 
interactions called from the background data; or using the assigned cutoffs; 4) 
estimate the significance based on the Poisson test for transformed data, both 
for the loop and loop anchors. For example, if transformed PETs for target is 
5, PETs for control is 3 while control nearby permutated background median is 
4, then for the Poisson test, lambda=4-1 is used to test the observed 5 to call
p-value.

Example:
    1. classical usage 
        cLoops2 callDiffLoops -tloop target_loop.txt -cloop control_loop.txt \\
                          -td ./target -cd ./control -o target_diff

    2. customize MA cutoffs 
        cLoops2 callDiffLoops -tloop target_loop.txt -cloop control_loop.txt \\
                          -td ./target -cd ./control -o target_diff -cutomize \\
                          -acut 5 -mcut 0.5
    """
    callDiffLoops = subparsers.add_parser(
        'callDiffLoops',
        parents=[parser],
        add_help=False,
        usage=argparse.SUPPRESS,
        description=callDiffLoopsDes,
        formatter_class=RawTextHelpFormatter,
    )
    callDiffLoops.add_argument(
        "-tloop",
        dest="tloop",
        type=str,
        required=True,
        help="The target loops in _loop.txt file called by cLoops2.")
    callDiffLoops.add_argument(
        "-cloop",
        dest="cloop",
        type=str,
        required=True,
        help="The control loops in _loop.txt file called by cLoops2.")
    callDiffLoops.add_argument(
        "-td",
        dest="tpred",
        type=str,
        required=True,
        help="The data directory generated by cLoops2 for target data.")
    callDiffLoops.add_argument(
        "-cd",
        dest="cpred",
        type=str,
        required=True,
        help="The data directory generated by cLoops2 for control data.")
    callDiffLoops.add_argument(
        "-pcut",
        dest="pcut",
        default=1e-2,
        type=float,
        help=
        "Poisson p-value cutoff to determine significant differentially\n"\
        "enriched loops after Bonferroni correction , default is 1e-2."
    )
    callDiffLoops.add_argument(
        "-igp",
        dest="igp",
        default=False,
        action="store_true",
        help=
        "Ignore Poisson p-value cutoff and only using FDR to control MA plot\n"\
        "cutoffs."
    )
    callDiffLoops.add_argument(
        "-noPCorr",
        dest="noPCorr",
        default=False,
        action="store_true",
        help=
        "Do not performe Bonferroni correction of Poisson p-values. Will get\n"\
        "more loops. Default is always performing."
    )
    callDiffLoops.add_argument(
        "-fdr",
        dest="fdr",
        default=0.1,
        type=float,
        help=
        "FDR cutoff for estimating fold change (M) and average value (A) after\n"\
        "normalization with background data. Default is 0.1."
    )
    callDiffLoops.add_argument(
        "-j",
        dest="juicebox",
        required=False,
        action="store_true",
        help=
        "Whether to convert loops to 2D feature annotations to visualize in\n"\
        "Juicebox. Default is not, set this flag to save."
    )
    callDiffLoops.add_argument(
        "-w",
        dest="washU",
        required=False,
        action="store_true",
        help=
        "Whether to save tracks of loops to visualize in legacy and new washU.\n"\
        "Default is not, set this flag to save two files."
    )
    callDiffLoops.add_argument(
        "-customize",
        dest="customize",
        required=False,
        action="store_true",
        help=
        "Whether to use cutomized cutoffs of MA plot. Defulat is not. If enable\n"\
        "-acut and -mcut is needed."
    )
    callDiffLoops.add_argument(
        "-cacut",
        dest="cacut",
        required=False,
        type=float,
        default=0.0,
        help=
        "Average cutoff for MA plot of normalized PETs of loops. Assign when\n"\
        "-customize option used."
    )
    callDiffLoops.add_argument(
        "-cmcut",
        dest="cmcut",
        required=False,
        type=float,
        default=0.0,
        help=
        "Fold change cutoff for MA plot of normalized PETs of loops. Assign when\n"\
        "-customize option used."
    )
    callDiffLoops.add_argument(
        "-vmin",
        dest="vmin",
        default=None,
        required=False,
        type=float,
        help=
        "The minimum value shown in the heatmap and colorbar."
    )
    callDiffLoops.add_argument(
        "-vmax",
        dest="vmax",
        default=None,
        required=False,
        type=float,
        help=
        "The maxmum value shown in the heatmap and colorbar."
    )
    callDiffLoops.add_argument(
        "-cmap",
        dest="cmap",
        default="summer",
        required=False,
        choices=["summer","red","div","cool"],
        help=
        "The heatmap style. Default is summer."
    )
 
    #call domain function
    callDomainsDes = """
Call domains for the 3D genomic data based on correlation matrix and local 
segregation score.

Well tested work for Hi-TrAC/Trac-looping data.

Examples:
    1. call Hi-C like TADs
        cLoops2 callDomains -d trac -o trac -bs 5000,10000 -ws 500000 -p 20

    2. call Hi-TrAC/Trac-looping specific small domains
        cLoops2 callDomains -d trac -o trac -bs 1000 -ws 100000 -p 20 

    3. call domains for Hi-C
        cLoops2 callDomains -d hic -o hic -bs 10000 -ws 500000 -hic 
"""
    callDomains = subparsers.add_parser(
        'callDomains',
        description=callDomainsDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )
    callDomains.add_argument(
        "-bs",
        dest="binSize",
        required=False,
        type=str,
        default="10000",
        help=
        "Candidate contact matrix resolution (bin size) to call domains. A\n"\
        "series of values or a single value can be used as input. Default is\n"\
        "10000. If given multiple values, callDomains will try to call nested\n"\
        "domains. Samll value may lead to samller domains."
    )
    callDomains.add_argument(
        "-ws",
        dest="winSize",
        required=False,
        type=str,
        default="500000",
        help=
        "The half of the sliding window size used to caculate local correlation,\n"\
        "Default is 500000 (500kb). If given multiple values, callDomains will\n"\
        "try to call nested domains. Larger value may lead to larger domains."
    )
    callDomains.add_argument(
        "-hic",
        dest="hic",
        required=False,
        action="store_true",
        help=
        "Whether to use cutoffs for Hi-C to output significant domains.\n"\
        "Default is not. Set this option to enable, cutoffs will be more loose."
    )
 
    #plot function
    plotDes = """
Plot the interaction data as a heatmap (or arches/scatter) with additional of 
virtual 4C view point, 1D tracks (bigWig files), 1D annotations (peaks, genes) 
and 2D annotations (domains). If -f is not assigned, will just plot profiles 
from bigWig file or bed files.

Examples:
    1. plot the simple square heatmap for a specific region with 1kb resolution 
       with genes 
        cLoops2 plot -f test/chr21-chr21.ixy -o test -bs 1000 -start 34840000 \\
                     -end 34895000 -log -gtf test.gtf

    2. plot the upper triangle heatmap with domains such as TAD and CTCF bigWig
       track
        cLoops2 plot -f test/chr21-chr21.ixy -o test_domain -bs 10000 \\
                     -start 34600000 -end 35500000 -domains HiC_TAD.bed -log \\
                    -triu -bws GM12878_CTCF_chr21.bw

    3. plot the heatmap as upper triangle with 1D signal track and filter the 
       PETs shorter than 1kb
        cLoops2 plot -f test/chr21-chr21.ixy -o test -bs 500 -start 34840000 \\
                     -end 34895000 -log -triu -1D -cut 1000

    4. plot the observation/expectation interaction heatmap with 1D signal 
        cLoops2 plot -f test/chr21-chr21.ixy -o test -m obs/exp -1D -triu \\ 
                     -bs 500 -start 34840000 -end 34895000

    5. plot the chromosome-wide correlation heatmap 
        cLoops2 plot -f test/chr21-chr21.ixy -o test -corr 

    6. plot upper triangle interaction heatmap together with genes, bigWig 
       files, peaks, loops, domains, control the heatmap scale
        cLoops2 plot -f test/chr21-chr21.ixy -o test -bs 500 -start 34840000 \\
                     -end 34895000 -triu -bws ATAC.bw,CTCF.bw -1D \\
                     -loop test_loops.txt -beds Enh.bed,Tss.bed \\
                     -domains tad.bed -m obs -log -vmin 0.2 -vmax 2 -gtf genes.gtf
    
    7. plot small regions interacting PETs as arches 
        cLoops2 plot -f test/chr21-chr21.ixy -o test -start 46228500 \\
                     -end 46290000 -1D -loops gm_loops.txt -arch -aw 0.05

    8. plot small regions interacting PETs as scatter plot
        cLoops2 plot -f test/chr21-chr21.ixy -o test -start 46228500 \\
                     -end 46290000 -1D -loops gm_loops.txt -scatter

    9. plot Hi-C compartments and eigenvector  
        cLoops2 plot -f test/chr21-chr21.ixy -o test -bs 100000 -log -corr -eig  
"""
    plot = subparsers.add_parser(
        'plot',
        description=plotDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )
    plot.add_argument(
        "-f",
        dest="fixy",
        required=False,
        type=str,
        default="",
        help=
        "Input .ixy file generated by cLoops2 pre. If not assigned, no heatmaps\n"\
        "or arches will be shown and -chrom is needed to generate plots similar\n"\
        "to IGV or other browser."
    )
    plot.add_argument(
        "-bs",
        dest="binSize",
        required=False,
        default=5000,
        type=int,
        help=
        "Bin size/matrix resolution (bp) to generate the contact matrix for\n"\
        "plotting, default is 5000 bp."
    )
    plot.add_argument(
        "-chrom",
        dest="chrom",
        required=False,
        type=str,
        default="",
        help="Chromosome for the target region if -f is not assigned."
    )
    plot.add_argument(
        "-start",
        dest="start",
        required=False,
        type=int,
        default=0,
        help="Start genomic coordinate for the target region. Default is 0."
    )
    plot.add_argument(
        "-end",
        dest="end",
        required=False,
        type=int,
        default=-1,
        help=
        "End genomic coordinate for the target region. Default is to infer\n"\
        "from the data."
    )
    plot.add_argument(
        "-loops",
        dest="floop",
        required=False,
        default="",
        type=str,
        help=
        "The _loop.txt file generated by cLoops2, will be used to plot loops as\n"\
        "arches."
    )
    plot.add_argument(
        "-domains",
        dest="fdomain",
        required=False,
        default="",
        type=str,
        help=
        "The domains need to annotated in the heatmap such as TADs, should be\n"\
        ".bed file."
    )
    plot.add_argument(
        "-beds",
        dest="beds",
        default="",
        required=False,
        help=
        "BED tracks of genomic features to plot above the heatmap, such as\n"\
        "promoters and enhancers, track name will be inferred from file name,\n"\
        "for example enhancer.bed,promoter.bed."
    )
    plot.add_argument(
        "-gtf",
        dest="gtf",
        default="",
        required=False,
        help=
        "GTF track of genes to plot above the heatmap."
    )
    plot.add_argument(
        "-bws",
        dest="bws",
        default="",
        required=False,
        help=
        "BigWig tracks to plot above the heatmap, track name will be inferred\n"\
        "from file name, for example a.bw,b.bw,c.bw. "
    )
    plot.add_argument(
        "-bwvs",
        dest="bwvs",
        default="",
        required=False,
        help=
        "BigWig tracks y-axis limitations. Default is atuo-determined. Assign\n"\
        "as 'vmin,vmax;vmin,vmax;vmin,vmax'. For example, '0,1;;0,1' for three\n"\
        "bigWig tracks, as the second track kept atuo-determined. Due to\n"\
        "argparse limitation for parsing minus value, also can be assigned as\n"\
        "vmax,vmin."
    )
    plot.add_argument(
        "-bwcs",
        dest="bwcs",
        default="",
        required=False,
        help=
        "BigWig tracks colors. Default is atuo-determined. Assign as \n"\
        "0,1,2 for three bigWig tracks. Values seperated by comma."
    )
    plot.add_argument(
        "-log",
        dest="log",
        required=False,
        action="store_true",
        default=False,
        help="Whether to log transform the matrix."
    )
    plot.add_argument(
        "-m",
        dest="method",
        type=str,
        choices=["obs", "obs/exp"],
        default="obs",
        help=
        "The type of matrix to plot, observed matrix or observed/expected\n"\
        "matrix, expected matrix will be generated by shuffling PETs, default\n"\
        "is observed."
    )
    plot.add_argument(
        "-corr",
        dest="corr",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to plot the correlation matrix. Default is not. Correlation\n"\
        "heatmap will use dark mode color map, used together with obs method."
    )
    plot.add_argument(
        "-norm",
        dest="norm",
        default=False,
        required=False,
        action="store_true",
        help="Whether to normalize the matrix with z-score."
    )
    plot.add_argument(
        "-triu",
        dest="triu",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to rotate the heatmap only show upper triangle, default is\n"\
        "False."
    )
    plot.add_argument(
        "-1D",
        dest="oneD",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to plot the pileup 1D signal for the region. Default is not.\n"\
        "Please note, the 1D signal is aggregated from the visualization region.\n"
        "If want to check the signal from each position of all genome/chromosome,\n"
        "use cLoops2 dump -bdg to get the bigWig file."
    )
    plot.add_argument(
        "-1Dv",
        dest="oneDv",
        default="",
        required=False,
        help=
        "1D profile y-axis limitations. Default is auto-determined. Assign as\n"\
        "vmin,vmax, for example 0,1."
    )
    plot.add_argument(
        "-vmin",
        dest="vmin",
        default=None,
        required=False,
        type=float,
        help=
        "The minimum value shown in the heatmap and colorbar."
    )
    plot.add_argument(
        "-vmax",
        dest="vmax",
        default=None,
        required=False,
        type=float,
        help=
        "The maxmum value shown in the heatmap and colorbar."
    )
    plot.add_argument(
        "-virtual4C",
        dest="virtual4C",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to plot the virtual 4C view point 1D signal. Default is not.\n"\
        "If assigned, -view_start and -view_end are needed."
    )
    plot.add_argument(
        "-view_start",
        dest="viewStart",
        required=False,
        type=int,
        default=0,
        help=
        "Start genomic coordinate for the view point start region, only valid\n"\
        "when -vitrutal4C is set, should >=start and <=end."
    )
    plot.add_argument(
        "-view_end",
        dest="viewEnd",
        required=False,
        type=int,
        default=-1,
        help=
        "End genomic coordinate for the view point end region, only valid\n"\
        "when -vitrutal4C is set, should >=start and <=end."
    )
    plot.add_argument(
        "-arch",
        dest="arch",
        required=False,
        action="store_true",
        default=False,
        help="Whether to plot interacting PETs as arches. Default is not. If\n"\
        "set, only original one PET one arch will be shown. Usefule to check\n"\
        "small region for raw data, especially when heatmap is not clear."
    )
    plot.add_argument(
        "-aw",
        dest="aw",
        required=False,
        default=1,
        type=float,
        help="Line width for each PET in arches plot. Default is 1. Try to\n"\
        "change it if too many or few PETs."
    )
    plot.add_argument(
        "-ac",
        dest="ac",
        required=False,
        default=4,
        type=int,
        help="Line color for each PET in arches plot. Default is 4. Try to\n"\
        "change it see how many colors are supported by cLoops2."
    )
    plot.add_argument(
        "-aa",
        dest="aa",
        required=False,
        default=1,
        type=float,
        help="Alpha to control arch color saturation. Default is 1."
    )
    plot.add_argument(
        "-scatter",
        dest="scatter",
        required=False,
        action="store_true",
        default=False,
        help="Whether to plot interacting PETs as scatter dots. Default is not.\n"\
        "If set, only original one PET one dot will be shown. Usefule to check\n"\
        "raw data, especially when heatmap is not clear that -vmax is too small."
    )
    plot.add_argument(
        "-ss",
        dest="ss",
        required=False,
        default=1,
        type=float,
        help="Dot size for each PET in scatter plot. Default is 1. Try to\n"\
        "change it to optimize the plot."
    )
    plot.add_argument(
        "-sc",
        dest="sc",
        required=False,
        default=0,
        type=int,
        help="Dot color for each PET in scatter plot. Default is 0. Try to\n"\
        "change it see how many colors are supported by cLoops2."
    )
    plot.add_argument(
        "-sa",
        dest="sa",
        required=False,
        default=1,
        type=float,
        help="Alpha to control dot color saturation. Default is 1."
    )
    plot.add_argument(
        "-eig",
        dest="eig",
        required=False,
        action="store_true",
        default=False,
        help="Whether to plot the PC1 of correlation matirx to show compartments\n"\
             "Default is not. Only work well for big regions such as resolution\n"\
             "of 100k."
    )
    plot.add_argument(
        "-eig_r",
        dest="eig_r",
        required=False,
        action="store_true",
        default=False,
        help="Whether to flip the PC1 values of -eig. It should be dependend on\n"\
        "inactivate or activate histone markers, as actually the PCA values do\n"\
        "not have directions, especially comparing different samples."
    )
    plot.add_argument(
        "-figWidth",
        dest="figWidth",
        required=False,
        default=4,
        choices=[4, 8],
        type=int,
        help="Figure width. 4 is good to show the plot as half of a A4 figure\n"\
        "width and 8 is good to show more wider. Default is 4."
    )

    #montage analysis
    montageDes = """
Montage analysis of specific regions, producing Westworld Season 3 -like 
Rehoboam plot. 

Examples: 
    1. showing all PETs for a gene's promoter and enhancers
        cLoops2 montage -f test/chr21-chr21.ixy -bed test.bed -o test 

    2. showing simplified PETs for a gene's promoter and enhancers
        cLoops2 montage -f test/chr21-chr21.ixy -bed test.bed -o test -simple
    
    3. adjust interacting link width 
        cLoops2 montage -f test/chr21-chr21.ixy -bed test.bed -o test -simple \\
                        -ppmw 10
    
    4. showing all PETs for a region, if in the bed file only contains one region
        cLoops2 montage -f test/chr21-chr21.ixy -bed test.bed -o test -ext 0
    """
    montage = subparsers.add_parser(
        'montage',
        description=montageDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )
    montage.add_argument(
        "-f",
        dest="fixy",
        required=True,
        type=str,
        help="Input .ixy file generated by cLoops2 pre."
    )
    montage.add_argument(
        "-bed",
        dest="bed",
        required=True,
        type=str,
        help=
        "Input .bed file for target regions, 4th columns should be id/name for\n"
        "the region."
    )
    montage.add_argument(
        "-ext",
        dest="ext",
        required=False,
        type=float,
        default=2,
        help=
        "Up-stream and down-stream extesion of target region length. Default is\n"\
        "2. If the input bed already include up/down-stream regions, assign as 0."
    )
    montage.add_argument(
        "-simple",
        dest="simple",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to only draw the representative interactions between two target\n"\
        "regions as one arch, and not include the interactions in extended\n"\
        "regions. Default is not, all interactions will be shown as archs.."
    )
    montage.add_argument(
        "-vp",
        dest="viewPoint",
        default="",
        type=str,
        required=False,
        help=
        "Only show interactions with specific regions from all other regions.\n"\
        "Name/id (4th column in .bed file) is need. Default is to show all\n"\
        "releated interactions. Multiple names/ids can be assigned by seperation\n"\
        "of comma."
    )
    montage.add_argument(
        "-vmin",
        dest="vmin",
        default=None,
        required=False,
        type=float,
        help=
        "The minial scale for 1D pileup data. Default will be inferred from the\n"\
        "data."
    )
    montage.add_argument(
        "-vmax",
        dest="vmax",
        default=None,
        required=False,
        type=float,
        help=
        "The maxmial scale for 1D pileup data. Default will be inferred from the\n"\
        "data."
    )
    montage.add_argument(
        "-ppmw",
        dest="ppmw",
        default=10,
        required=False,
        type=float,
        help=
        "Link line width indicator, short for 1 PETs per Million PETs line\n"
        "width, default is 10. Adjust this value when -simple is used. Decrease\n"\
        "it if links are too bold and increase it when links are too thin."
    )
    montage.add_argument(
        "-aw",
        dest="aw",
        required=False,
        default=1,
        type=float,
        help="Line width for each PET if -simple is not selected. Default is 1."
    )
    montage.add_argument(
        "-no1D",
        dest="noOneD",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to not plot 1D profiles. Default is plot. Set this for Hi-C\n"\
        "like data."
    )
 
 
    aggDes = """
Do the aggregation analysis for peaks, loops, view points and domains.

The output figures can be used directly, and the data to generate the plot are 
also saved for further customized analysis. 

For the aggregated peaks analysis,input is a .bed file annotated with the 
coordinates for the target regions/peaks/anchors. Output is a .pdf file 
containing a mean density plot and heatmap and a .txt file for the data. The 
data in the .txt file and plot were normalized to RPM (reads per million).

For the aggregated view points analysis, input is a .bed file annotated with 
coordinates for the target regions/peaks/anchors as view point. Output is a 
.pdf file containing a mean density plot and heatmap and a .txt file for the 
data. The data in the .txt file and plot were normalized to 
log2( RPM (reads per million)+1).

For the aggregated loops analysis, input is a _loops.txt file annotated with 
the coordinates for target loops, similar to the format of BEDPE. Output is a 
.pdf file for mean heatmap and .npz file generated through numpy.savez for all 
loops and nearby regions matrix. The enrichment score (ES) in the plot is 
calculated as: ES = mean( (PETs in loop)/(mean PETs of nearby regions) ). Other 
files except _loops.txt can be used as input, as long as the file contains key 
information in the first columns separated by tabs:
loopId\tchrA\tstartA\tendA\tchrB\tstartB\tendB\tdistance
loop-1\tchr21\t1000\t2000\tchr21\t8000\t9000\t7000\n

There is another option for loops analysis, termed as two anchors. Input file is 
same to aggregated loops analysis. The whole region with assigned extesion
between two anchors will be aggregated and 1D profile can show two anchors. The 
analysis could be usefule to study/comapre different classes of anchors and 
combinations, for example, considering CTCT motif directions, all left anchors
CTCF motifs are in positive strand and in negative strand for all right anchors. 
It could be interesting for some loops one anchor only bound by transcription 
factor a and another anchor only bound by transcription b. 

For the aggregated domains analysis, input is a .bed file annotated with the
coordinates for the domains, such as TADs. Output are a .pdf file for the upper 
triangular heatmap and .npz file generated through numpy.savez for all domains 
and nearby region matrix. The enrichment score (ES) in the plot is calculated 
as mean( (two ends both with in domain PETs number)/( only one end in domain 
PETs number) ).

Examples:
    1. show aggregated peaks heatmap and profile 
        cLoops2 agg -d test -peaks peaks.bed -o test -peak_ext 2500 \\ 
                    -peak_bins 200 -peak_norm -skipZeros

    2. show aggregated view points and aggregated bigWig signal
        cLoops2 agg -d test -o test -viewPoints test_peaks.bed -bws CTCF.bw 

    3. show aggregated loops heatmap, 1D profile and aggregated bigWig signal
        cLoops2 agg -d test -o test -loops test_loops.txt -bws CTCF.bw -1D \\
                    -loop_norm
    
    3. show aggregated loops heatmap, 1D profile and aggregated bigWig signal
       in two anchors mode
        cLoops2 agg -d test -o test -twoAnchors test_loops.txt -bws CTCF.bw -1D \\
                    -loop_norm

    4. show aggregated domains heatmap, 1D profile and aggregated bigWig signal
        cLoops2 agg -d test -o test -domains TAD.bed -bws CTCF.bw -1D 
    """
    agg = subparsers.add_parser(
        'agg',
        description=aggDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )
    agg.add_argument(
        "-peaks",
        dest="peakf",
        required=False,
        default="",
        type=str,
        help="The .bed file for peaks-centric aggregation analysis."
    )
    agg.add_argument(
        "-peak_ext",
        dest="peak_ext",
        required=False,
        default=5000,
        type=int,
        help=
        "The nearby upstream and downstream regions (bp) from the peak center.\n"\
        "Default is 5000."
    )
    agg.add_argument(
        "-peak_bins",
        dest="peak_bins",
        required=False,
        default=100,
        type=int,
        help="The bin size for the profile array of peaks. Default is 100."
    )
    agg.add_argument(
        "-peak_norm",
        dest="peak_norm",
        required=False,
        default=False,
        action="store_true",
        help=
        "Whether to normalize the data in the peaks profile plot and\n"\
        "heatmap with row-wise z-score. Default is not."
    )
    agg.add_argument(
        "-viewPoints",
        dest="viewPointF",
        required=False,
        default="",
        type=str,
        help="The .bed file for view points -centric aggregation analysis."
    )
    agg.add_argument(
        "-viewPointUp",
        dest="viewPointUp",
        required=False,
        default=100000,
        type=int,
        help=
        "The upstream regions included for the aggreaged view points analysis.\n"\
        "Default is 100000 bp."
    )
    agg.add_argument(
        "-viewPointDown",
        dest="viewPointDown",
        required=False,
        default=100000,
        type=int,
        help=
        "The downstream regions included for the aggreaged view points analysis.\n"\
        "Default is 100000 bp."
    )
    agg.add_argument(
        "-viewPointBs",
        dest="viewPointBs",
        required=False,
        default=1000,
        type=int,
        help=
        "Contact matrix bin size for view points heatmap. Default is 1000 bp. "
    )
    agg.add_argument(
        "-viewPoint_norm",
        dest="viewPoint_norm",
        required=False,
        default=False,
        action="store_true",
        help=
        "Whether to normalize the sub-matrix for each loop as divide the mean\n"\
        "PETs for the matrix. Default is not."
    )
    agg.add_argument(
        "-viewPoint_vmin",
        dest="viewPoint_vmin",
        required=False,
        default=None,
        type=float,
        help=
        "The minimum value shown in the aggregated view points heatmap and colorbar."
    )
    agg.add_argument(
        "-viewPoint_vmax",
        dest="viewPoint_vmax",
        required=False,
        default=None,
        type=float,
        help=
        "The maxmum value shown in the aggregated view points heatmap and colorbar."
    )
    agg.add_argument(
        "-loops",
        dest="loopf",
        required=False,
        default="",
        type=str,
        help=
        "The _loop.txt file generated by cLoops2 for loops-centric\n"\
        "aggregation analysis. The file first 8 columns are necessary."
    )
    agg.add_argument(
        "-loop_ext",
        dest="loop_ext",
        required=False,
        default=10,
        type=int,
        help=
        "The nearby regions included to plot in the heatmap and calculation of\n"
        "enrichment for aggregation loop analysis, default is 10, should be\n"\
        "even number."
    )
    agg.add_argument(
        "-loop_cut",
        dest="loop_cut",
        type=int,
        default=0,
        help="Distance cutoff for loops to filter. Default is 0."
    )
    agg.add_argument(
        "-loop_norm",
        dest="loop_norm",
        required=False,
        default=False,
        action="store_true",
        help=
        "Whether to normalize the sub-matrix for each loop as divide the mean\n"\
        "PETs for the matrix (except the loop region). Default is not."
    )
    agg.add_argument(
        "-loop_vmin",
        dest="loop_vmin",
        required=False,
        default=None,
        type=float,
        help=
        "The minimum value shown in the aggregated loops heatmap and colorbar."
    )
    agg.add_argument(
        "-loop_vmax",
        dest="loop_vmax",
        required=False,
        default=None,
        type=float,
        help=
        "The maxmum value shown in the aggregated loops heatmap and colorbar."
    )
    agg.add_argument(
        "-twoAnchors",
        dest="twoAnchorsF",
        required=False,
        default="",
        type=str,
        help=
        "The similar _loop.txt file generated by cLoops2 for two anchors\n"\
        "aggregation analysis. The file first 8 columns are necessary."
    )
    agg.add_argument(
        "-twoAnchor_ext",
        dest="twoAnchor_ext",
        required=False,
        default=0.1,
        type=float,
        help="The nearby regions of fold included to plot in heatmap.\n"\
        "Default is 0.1."
    )
    agg.add_argument(
        "-twoAnchor_vmin",
        dest="twoAnchor_vmin",
        default=None,
        required=False,
        type=float,
        help=
        "The minimum value shown in the domain heatmap and colorbar."
    )
    agg.add_argument(
        "-twoAnchor_vmax",
        dest="twoAnchor_vmax",
        default=None,
        required=False,
        type=float,
        help=
        "The maxmum value shown in the domain heatmap and colorbar."
    )
    agg.add_argument(
        "-domains",
        dest="domainf",
        required=False,
        default="",
        type=str,
        help="The .bed file annotated the domains such as TADs for aggregated\n"\
        "domains-centric analysis."
    )
    agg.add_argument(
        "-domain_ext",
        dest="domain_ext",
        required=False,
        default=0.5,
        type=float,
        help="The nearby regions of fold included to plot in heatmap and\n"\
        "caculation of enrichment, default is 0.5."
    )
    agg.add_argument(
        "-domain_vmin",
        dest="domain_vmin",
        default=None,
        required=False,
        type=float,
        help=
        "The minimum value shown in the domain heatmap and colorbar."
    )
    agg.add_argument(
        "-domain_vmax",
        dest="domain_vmax",
        default=None,
        required=False,
        type=float,
        help=
        "The maxmum value shown in the domain heatmap and colorbar."
    )
    agg.add_argument(
        "-1D",
        dest="oneD",
        default=False,
        required=False,
        action="store_true",
        help="Whether to plot the pileup 1D signal for aggregated loops, \n"\
        "aggregated view points or aggregated domains. Default is not."
    )
    agg.add_argument(
        "-bws",
        dest="bws",
        default="",
        required=False,
        help=
        "BigWig tracks to plot above the aggregated loops heatmap (or under\n"\
        "the aggregated domains heatmap), track name will be inferred from file\n"\
        "name, for example a.bw,b.bw,c.bw. "
    )
    agg.add_argument(
        "-skipZeros",
        dest="skipZeros",
        required=False,
        default=False,
        action="store_true",
        help="Whether to remove all 0 records. Default is not."
    )

    quantDes = """
Quantify the peaks, loops and domains.  The output file will be the same as
outputs of callPeaks, callLoops and callDomains.

Examples:
    1. quantify peaks 
        cLoops2 quant -d test -peaks peaks.bed -o test 

    2. quantify loops 
        cLoops2 quant -d test -loops test_loops.txt -o test
    
    3. quantify domains 
        cLoops2 quant -d test -domains test_domains.txt -o test
"""
    quant = subparsers.add_parser(
        'quant',
        description=quantDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )
    quant.add_argument(
        "-peaks",
        dest="peakf",
        required=False,
        default="",
        type=str,
        help="The .bed file for peaks-centric quantification."
    )
    quant.add_argument(
        "-loops",
        dest="loopf",
        required=False,
        default="",
        type=str,
        help="The _loop.txt file generated by cLoops2 for loops-centric\n"\
        "quantification, as long as there are first 8 columns."
    )
    quant.add_argument(
        "-domains",
        dest="domainf",
        required=False,
        default="",
        type=str,
        help="The _domains.txt file generated by cLoops2 for domains-centric\n"\
        "quantification, as long as there are first 3 columns"
    )
    quant.add_argument(
        "-domain_bs",
        dest="domainBinSize",
        required=False,
        type=int,
        default=10000,
        help=
        "Candidate contact matrix resolution (bin size) to quantify domains, \n"\
        "default is 10000. Only one integer is supported.\n"
    )
    quant.add_argument(
        "-domain_ws",
        dest="domainWinSize",
        required=False,
        type=int,
        default=500000,
        help=
        "The half window size used to calculate local correlation to quantify\n"\
        "domains. Default is 500000 (500kb)."
    )
    quant.add_argument(
        "-domain_bdg",
        dest="domainBdg",
        default=False,
        required=False,
        action="store_true",
        help="Whether to save the segregation score as bedGraph file, default\n"\
        "is not."
    )

    anaLoopsDes = """
Annotating loops:
- find the closest TSS for each loop anchors
- merge the loop anchors and classify them as enhancers or promoters based on 
  distance to nearest TSS
- build the interaction networks for merged anchors 
- find the all interacted enhancers/promoters for each promoter  

Basic mode 1: with -gtf, loops will be annotated as enhancer or promoter based 
on distance to nearest gene. If a anchor overlapped with two/multiple promoters
(often seen for close head-to-head genes), all will be reported. If no TSS 
overlaps, then nearest one will be assigned.  

Basic mode 2: with -gtf -net, overlapped anchors will be merged and annoated as 
enhancer or promoter considering distance to genes. For each promoter, all 
linked enhancer and promoter will be shown. If there are more than 3 direct or 
indirect enhancers for a promoter, HITS algorithm will be used to identify one
hub for indirect enhancer and one hub for indirect enhancer. 

Examples:
    1. annotate loops for target gene, basic mode 1
        cLoops2 anaLoops -loops test_loops.txt -gtf genecode.gtf
    
    2. annotate loops for target transcripts (alternative TSS), basic mode 1
        cLoops2 anaLoops -loops test_loops.txt -gtf genecode.gtf -tid
    
    3. find a gene's all linked enhancer or promoter, basic mode 2
        cLoops2 anaLoops -loops test_loops.txt -gtf genecode.gtf -net
"""
    anaLoops = subparsers.add_parser(
        'anaLoops',
        description=anaLoopsDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )    
    anaLoops.add_argument(
        "-loops",
        dest="floop",
        required=True,
        default="",
        type=str,
        help=
        "The _loop.txt file generated by cLoops2 callLoops or callDiffLoops."
    )
    anaLoops.add_argument(
        "-gtf",
        dest="gtf",
        default="",
        required=False,
        type=str,
        help=
        "GTF file annotation for genes."
    )
    anaLoops.add_argument(
        "-tid",
        dest="tid",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to use transcript id instead of gene id for annotation. Default\n"\
        "is not."
    )
    anaLoops.add_argument(
        "-pdis",
        dest="pdis",
        default=2000,
        required=False,
        type=int,
        help=
        "Distance limitation for anchor to nearest gene/transcript TSS to define\n"\
        "as promoter. Default is 2000 bp."
    )
    anaLoops.add_argument(
        "-net",
        dest="net",
        default=False,
        required=False,
        action="store_true",
        help=
        "Whether to use network method to find all enhancer/promoter links based\n"\
        "on loops. Default is not. In this mode, overlapped anchors will be\n"\
        "merged and annotated as enhancer/promoter, then for a gene, all linked\n"\
        "node will be output."
    )
    anaLoops.add_argument(
        "-gap",
        dest="gap",
        default=1,
        required=False,
        type=int,
        help=
        "When -net is set, the distance for close anchors to merge. Default is 1."
    )
   
    findTargetsDes = """
Find target genes of genomic regions (peaks, SNPs) through enhancer-promoter 
networks. Output from cLoops2 anaLoops with suffix of _ep_net.sif and
_targets.txt are needed.

Examples:
    1. find target genes of peaks/SNPs
        cLoops2 findTargets -net test_ep_net.sif -tg test_targets.txt \\
                            -bed GWAS.bed -o test 
"""
    findTargets = subparsers.add_parser(
        'findTargets',
        description=findTargetsDes,
        formatter_class=RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        parents=[parser],
        add_help=False,
    )    
    findTargets.add_argument(
        "-net",
        dest="fnet",
        required=True,
        default="",
        type=str,
        help=
        "The _ep_net.sif file generated by cLoops2 anaLoops."
    )
    findTargets.add_argument(
        "-tg",
        dest="ftg",
        required=True,
        default="",
        type=str,
        help=
        "The _targets.txt file generated by cLoops2 anaLoops."
    )
    findTargets.add_argument(
        "-bed",
        dest="fbed",
        required=False,
        default="",
        help= "Find target genes for regions, such as anchors, SNPs or peaks."
    )

    return parser

def main():
    cliParser = mainHelp()
    if len(sys.argv) < 2:  #no subcommand assigned, just show help  and return
        cliParser.print_help()
        return
    #get the specific command
    cmd = sys.argv[1]
    #get all parser
    cliParser = cliParser.parse_args()

    #1. quality control
    if cmd == "qc":
        start = datetime.now()

        report = "Command: cLoops2 {} -f {} -o {} -p {} ".format(
            cmd, cliParser.fnIn, cliParser.fnOut, cliParser.cpu)
        logger.info(report)

        fs = cliParser.fnIn.split(",")
        for f in fs:
            if not os.path.isfile(f):
                report = "Input file %s not exitst!" % f
                logger.error(report)
                return
        #check output files
        fout = cliParser.fnOut + "_bedpeQc.txt"
        if os.path.isfile(fout):
            logger.error("Output file %s exists, return.\n" % fout)
            return

        qcBedpes(fs, fout, cliParser.cpu)

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #2. pre-process PETs into cLoops2 data
    if cmd == "pre":
        start = datetime.now()

        report = "Command: cLoops2 {} -f {} -o {} -p {} -c {} -cut {} -mcut {} -mapq {} -trans {}".format(
                cmd, 
                cliParser.fnIn, 
                cliParser.fnOut, 
                cliParser.cpu,
                cliParser.chroms, 
                cliParser.cut, 
                cliParser.mcut,
                cliParser.mapq, 
                cliParser.trans
            )
        logger.info(report)

        fout = cliParser.fnOut
        if not os.path.isdir(fout):
            os.mkdir(fout)
        elif len(os.listdir(fout)) > 0:
            r = "working directory %s exists and not empty." % fout
            logger.error(r)
            return

        if cliParser.chroms == "":
            chroms = []
        else:
            chroms = set(cliParser.chroms.split(","))

        fs = cliParser.fnIn.split(",")
        for f in fs:
            if not os.path.isfile(f):
                report = "Input file %s not exitst!" % f
                logger.error(report)
                return
        #parse BEDPE files into xy coordinates
        if cliParser.trans:
            cis = False
        else:
            cis = True
        parseBedpe(fs,
                   fout,
                   logger,
                   mapq=cliParser.mapq,
                   cs=chroms,
                   cut=cliParser.cut,
                   mcut=cliParser.mcut,
                   cis=cis,
                   cpu=cliParser.cpu,
                   )
        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #3. update .ixy and .json files just in case the directory is moved or manuplicated
    if cmd == "update":
        start = datetime.now()

        report = "Command: cLoops2 {} -d {}".format(cmd, cliParser.predir)
        logger.info(report)

        #check the input directory
        if cliParser.predir == "":
            logger.error("-d is required, return.")
            return
        writeNewJson(cliParser.predir)

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)
    
    #4. combine mulitple cLoops2 datasets 
    if cmd == "combine":
        start = datetime.now()

        report = "Command: cLoops2 {} -ds {} -o {} -keep {} -cpu {}".format(
                    cmd, 
                    cliParser.predirs,
                    cliParser.fnOut, 
                    cliParser.keep,
                    cliParser.cpu,
                    )
        logger.info(report)

        #check the input directory
        if cliParser.predirs is None or cliParser.predirs == "":
            logger.error("-ds is required, return.")
            return
        ds = cliParser.predirs.split(",")
        if len(ds) < 2: 
            logger.error("More than one sample shoud be assigned, return.")
            return
        flag = False
        for d in ds:
            if not os.path.exists(d):
                logger.error("%s not exists."%d)
                flag = True
        if flag:
            logger.error("Above inputs not exists, return.")
            return 
        if os.path.exists( cliParser.fnOut):
            logger.error("Output directory %s exists, return."%cliParser.fnOut)
            return
         
        combineDirs(    
            ds,
            cliParser.fnOut,
            logger,
            keep=cliParser.keep,
            cpu=cliParser.cpu
        )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)
    
    #5. dump cLoops2 data into other kinds of format
    if cmd == "dump":
        start = datetime.now()

        report = "Command: cLoops2 {cmd} -d {predir} -o {output} -cut {cut} -mcut {mcut} -bed {bed} -bed_ext {bed_ext} -bedpe {bedpe} -bedpe_ext {bedpe_ext} -hic {hic} -hic_org {hic_org} -hic_res {hic_res} -washU {washu} -washU_ext {washu_ext} -ucsc {ucsc} -ucsc_ext {ucsc_ext} -ucsc_cs {ucsc_cs} -bdg {bdg} -bdg_ext {bdg_ext} -bdg_pe {bdg_pe} -mat {mat} -mat_res {mat_res} -mat_chrom {chrom} -mat_start {start} -mat_end {end} -log {log} -m {method} -corr {corr} -norm {norm}".format(
            cmd = cmd, 
            predir = cliParser.predir,
            output = cliParser.fnOut,
            cut = cliParser.cut,
            mcut = cliParser.mcut,
            bed = cliParser.bed,
            bed_ext = cliParser.bed_ext,
            bedpe = cliParser.bedpe,
            bedpe_ext = cliParser.bedpe_ext,
            hic = cliParser.hic,
            hic_org = cliParser.hic_org,
            hic_res = cliParser.hic_res,
            washu = cliParser.washU,
            washu_ext = cliParser.washU_ext,
            ucsc = cliParser.ucsc,
            ucsc_ext = cliParser.ucsc_ext,
            ucsc_cs = cliParser.ucsc_cs,
            bdg = cliParser.bdg,
            bdg_ext = cliParser.bdg_ext,
            bdg_pe = cliParser.bdg_pe,
            mat = cliParser.mat,
            mat_res = cliParser.mat_res,
            chrom = cliParser.chrom,
            start = cliParser.start,
            end = cliParser.end,
            log = cliParser.log,
            method = cliParser.method,
            corr = cliParser.corr,
            norm = cliParser.norm,
        )
        logger.info(report)

        #check the input directory
        if cliParser.predir == "":
            logger.error("-d is required, return.")
            return

        #convert to .bed 
        if cliParser.bed:
            ixy2bed(
                cliParser.predir,
                cliParser.fnOut,
                logger,
                cut=cliParser.cut,
                mcut=cliParser.mcut,
                ext=cliParser.bed_ext,
            )


        #convert to .bedpe 
        if cliParser.bedpe:
            ixy2bedpe(
                cliParser.predir,
                cliParser.fnOut,
                logger,
                cut=cliParser.cut,
                mcut=cliParser.mcut,
                ext=cliParser.bedpe_ext,
            )

        #convert to .hic 
        if cliParser.hic:
            ixy2hic(
                cliParser.predir,
                cliParser.fnOut,
                logger,
                org=cliParser.hic_org,
                resolution=cliParser.hic_res,
                cut=cliParser.cut,
                mcut=cliParser.mcut,
            )
        #convert to washU track
        if cliParser.washU:
            ixy2washU(
                cliParser.predir,
                cliParser.fnOut,
                logger,
                cut=cliParser.cut,
                mcut=cliParser.mcut,
                ext=cliParser.washU_ext,
            )
        #convert to washU track
        if cliParser.ucsc:
            if not os.path.isfile( cliParser.ucsc_cs ):
                logger.error("A chrom sizes file needed, %s not exists, return."%cliParser.ucsc_cs)
                return
            ixy2ucsc(
                cliParser.predir,
                cliParser.fnOut,
                cliParser.ucsc_cs,
                logger,
                cut=cliParser.cut,
                mcut=cliParser.mcut,
                ext=cliParser.ucsc_ext,
            )
        #convert to bedGraph
        if cliParser.bdg:
            ixy2bdg(    
                cliParser.predir,
                cliParser.fnOut,
                logger,
                cut=cliParser.cut,
                mcut=cliParser.mcut,
                ext=cliParser.bdg_ext,
                pe=cliParser.bdg_pe,
            )
        #convert to matrix txt file
        if cliParser.mat:
            ixy2mat(
                cliParser.predir,
                cliParser.fnOut,
                logger,
                chrom=cliParser.chrom,
                start=cliParser.start,
                end=cliParser.end,
                r=cliParser.mat_res,
                cut=cliParser.cut,
                mcut=cliParser.mcut,
                log=cliParser.log,
                method=cliParser.method,
                corr=cliParser.corr,
                norm=cliParser.norm,
            )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #6. estimate eps
    if cmd == "estEps":
        start = datetime.now()

        report = "Command: cLoops2 {} -d {} -fixy {} -o {} -p {} -method {} -k {} -cut {} -mcut {}".format(
                cmd,
                cliParser.predir, 
                cliParser.fixy, 
                cliParser.fnOut,
                cliParser.cpu, 
                cliParser.epsMethod, 
                cliParser.knn, 
                cliParser.cut, 
                cliParser.mcut
            )
        logger.info(report)

        #check output file
        fout = cliParser.fnOut + "_" + cliParser.epsMethod + "_eps.pdf"
        if os.path.isfile(fout):
            r = "Output file %s exists, return." % fout
            logger.error(r)
            return

        #check the input directory and file
        fs = glob(os.path.join(cliParser.predir, "*.ixy"))
        f = cliParser.fixy
        if len(fs) == 0:
            if not os.path.isfile(cliParser.fixy):
                r = "No input directory or file!"
                logger.error(r)
                return
            else:
                #use the file to estimate eps
                fs = [f]
        else:
            if os.path.isfile(cliParser.fixy):
                r = "Both input directory and file exists! Just input one!"
                logger.error(r)
                return
            else:
                #use the files to estimate eps
                fs = fs

        #check the methods of estimation
        if cliParser.epsMethod == "gmm":
            #get all distance (abs(Y-X))
            dis = []
            ds = Parallel(n_jobs=cliParser.cpu,backend="multiprocessing")(
                delayed(getXyDis)(f, cliParser.cut, cliParser.mcut) for f in fs)
            for d in ds:
                if d is not None:
                    dis.extend(d)
            dis = np.log2(np.array(dis))
            #estimate kinds of points and eps
            ps, eps = getGmmLabelsEps(dis)
            #plot the estimation
            plotGmmEst(dis, ps, eps, fout)

        elif cliParser.epsMethod == "k-dis":
            #check if cliParser.knn assigned
            if cliParser.knn == 0:
                r = "-k not assigned!"
                logger.error(r)
                return

            #get the k-distance
            dis = []
            ds = Parallel(n_jobs=cliParser.cpu,backend="multiprocessing")(
                delayed(getKDis)(f, cliParser.knn, cliParser.cut,cliParser.mcut) for f in fs)
            for d in ds:
                if d is not None:
                    dis.extend(d)
            dis = np.log10(np.array(dis))
            #sort accending
            dis = np.sort(dis)
            knee, eps = getKDisKneeEps(dis)
            #plotKDisE( dis, cliParser.knn, knee,eps, fout)
            plotKDis(
                dis, cliParser.knn,
                fout)  #currently maybe only visual check by human is better
        else:
            logger.error("Method not implemented,return.")
            return

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #7. estimate resolution
    if cmd == "estRes":
        start = datetime.now()

        report = "Command: cLoops2 {} -d {} -o {} -p {} -bs {} -cut {} -mcut {}".format(
                cmd, 
                cliParser.predir, 
                cliParser.fnOut, 
                cliParser.cpu,
                cliParser.binSize, 
                cliParser.cut,
                cliParser.mcut
            )
        logger.info(report)

        #check the input directory
        if cliParser.predir == "":
            logger.error("-d is required, return.")
            return

        cliParser.binSize = parseEps(cliParser.binSize)

        #do the job
        estRes(cliParser.predir,
               cliParser.fnOut,
               logger,
               bs=cliParser.binSize,
               cpu=cliParser.cpu,
               cut=cliParser.cut,
               mcut=cliParser.mcut,
               )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #8. estimate interaction distance limitation
    if cmd == "estDis":
        start = datetime.now()

        report = "Command: cLoops2 {cmd} -d {d} -o {fout} -p {cpu} -bs {bs} -c {chrom} -cut {cut} -mcut {mcut} -r {r} -plot {plot}".format(
                cmd = cmd, 
                d = cliParser.predir, 
                fout = cliParser.fnOut, 
                cpu = cliParser.cpu,
                bs = cliParser.binSize, 
                chrom = cliParser.chroms,
                cut = cliParser.cut,
                mcut = cliParser.mcut,
                r = cliParser.repeats,
                plot = cliParser.plot
            )
        logger.info(report)
        #check the input directory
        if cliParser.predir == "":
            logger.error("-d is required, return.")
            return
        #do the job
        estDis(cliParser.predir,
               cliParser.fnOut,
               bs=cliParser.binSize,
               cpu=cliParser.cpu,
               cut=cliParser.cut,
               mcut=cliParser.mcut,
               chroms=cliParser.chroms,
               repeats=cliParser.repeats,
               plot=cliParser.plot
               )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)


    #9. estimate sequencing saturation 
    if cmd == "estSat":
        start = datetime.now()

        report = "Command: cLoops2 {} -d {} -o {} -p {} -bs {} -tol {} -cut {} -mcut {}".format(
                cmd, 
                cliParser.predir, 
                cliParser.fnOut, 
                cliParser.cpu,
                cliParser.binSize, 
                cliParser.tol,
                cliParser.cut,
                cliParser.mcut
            )
        logger.info(report)

        #check the input directory
        if cliParser.predir == "":
            logger.error("-d is required, return.")
            return

        cliParser.binSize = parseEps(cliParser.binSize)
        #do the job
        estSat(cliParser.predir,
               cliParser.fnOut,
               logger,
               bs=cliParser.binSize,
               tol=cliParser.tol,
               cpu=cliParser.cpu,
               cut=cliParser.cut,
               mcut=cliParser.mcut,
               )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #10. estimate similarity
    if cmd == "estSim":
        start = datetime.now()

        report = "Command: cLoops2 {} -ds {} -o {} -cpu {} -cut {} -mcut {} -bs {} -m {} -pcut {} -n {} -plot {}".format(
             cmd, 
             cliParser.predirs,
             cliParser.fnOut, 
             cliParser.cpu,
             cliParser.cut,
             cliParser.mcut,
             cliParser.binSize,
             cliParser.method,
             cliParser.pcut,
             cliParser.n_components,
             cliParser.plot,
        )

        logger.info(report)
        #check the input directory
        if cliParser.predirs is None or cliParser.predirs == "":
            logger.error("-ds is required, return.")
            return
        ds = cliParser.predirs.split(",")
        if len(ds) < 2: 
            logger.error("More than one sample shoud be assigned, return.")
            return
        flag = False
        for d in ds:
            if not os.path.exists(d):
                logger.error("%s not exists."%d)
                flag = True
        if flag:
            logger.error("Above inputs not exists, return.")
            return 
 
        #do the job
        estSim(
            ds,
            cliParser.fnOut,
            bs=cliParser.binSize,
            method=cliParser.method,
            cut=cliParser.cut,
            mcut=cliParser.mcut,
            cpu=cliParser.cpu,
            pcut=cliParser.pcut,
            n_comp=cliParser.n_components,
            plot=cliParser.plot,
        )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #11. filter PETs
    if cmd == "filterPETs":
        start = datetime.now()
        report = "Command: cLoops2 {} -d {} -peak {} -loop {} -p {} -o {} -gap {} -singleton {} -bs {} -knn {} -eps {} -minPts {} -iv {} -both {}".format(
            cmd, cliParser.predir, cliParser.fbed, cliParser.floop,
            cliParser.cpu, cliParser.fnOut, cliParser.gap, cliParser.singleton, 
            cliParser.binSize, cliParser.knn, cliParser.eps, 
            cliParser.minPts, cliParser.iv, cliParser.both)
        logger.info(report)

        #check the input directory
        if cliParser.predir == "":
            logger.error("-d is required, return.")
            return

        #logic judgement
        ts = [ ]
        if cliParser.fbed != "" and os.path.isfile(cliParser.fbed):
            ts.append( 1 )
        else:
            ts.append( 0 )
        if cliParser.floop != "" and os.path.isfile(cliParser.floop):
            ts.append( 1 )
        else:
            ts.append( 0 )
        if cliParser.singleton:
            ts.append( 1 )
        else:
            ts.append( 0 )
        if cliParser.knn:
            ts.append( 1 )
        else:
            ts.append( 0 )

        if sum(ts) != 1:
            r = "No filtering option selected or multiple filtering options selected, only one is allowed. Return."
            logger.error(r)
            return

        #prepare output directory
        foutdir = cliParser.fnOut
        if not os.path.exists(foutdir):
            os.mkdir(foutdir)
        elif len(os.listdir(foutdir)) > 0:
            r = "Working directory %s exists and not empty. Return." % foutdir
            logger.error(r)
            return

        #filter by peaks
        if cliParser.fbed != "" and os.path.isfile(cliParser.fbed):
            filterPETsByPeaks(cliParser.predir, cliParser.fbed, foutdir,
                              cliParser.cpu, cliParser.iv, cliParser.gap)

        #filter by loops
        if cliParser.floop != "" and os.path.isfile(cliParser.floop):
            filterPETsByLoops(cliParser.predir, cliParser.floop, foutdir,
                              cliParser.cpu, iv=cliParser.iv, gap=cliParser.gap,
                              both=cliParser.both
                              )

        #filter by singleton
        if cliParser.singleton == True:
            if cliParser.binSize <= 0:
                logger.error(
                    "-bs assigned as %s, not valid, should be > 0 return." %
                    cliParser.binSize)
            filterPETsBySingletons( cliParser.predir, foutdir, 
                                    cliParser.binSize, cliParser.cpu)
        
        #filter by KNNs
        if cliParser.knn == True:
            filterPETsByKNNs( cliParser.predir, foutdir, 
                            eps=cliParser.eps, minPts=cliParser.minPts, 
                            cpu=cliParser.cpu )

        end = datetime.now()
        logger.info("cLoops2 %s finished for %s. Used time: %s." %
                    (cmd, cliParser.fnOut, end - start) + "\n" * 3)

    #12. sample PETs
    if cmd == "samplePETs":
        start = datetime.now()
        report = "Command: cLoops2 {} -d {} -o {} -tot {} -p {} ".format(
            cmd, cliParser.predir, cliParser.fnOut, cliParser.tot,
            cliParser.cpu)
        logger.info(report)

        #check the input directory
        if cliParser.predir == "":
            logger.error("-d is required, return.")
            return

        #prepare output directory
        foutdir = cliParser.fnOut
        if not os.path.exists(foutdir):
            os.mkdir(foutdir)
        elif len(os.listdir(foutdir)) > 0:
            r = "Working directory %s exists and not empty. Return." % foutdir
            logger.error(r)
            return

        #do the job
        samplePETs(cliParser.predir, foutdir, cliParser.tot, cpu=cliParser.cpu)

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #13. call peaks
    if cmd == "callPeaks":
        start = datetime.now()
        report = "Command: cLoops2 {cmd} -d {pdir} -cut {cut} -mcut {mcut} -eps {eps} -minPts {minPts} -pcut {pcut} -p {cpu} -o {output} -bgd {bgd} -bgm {bgm} -pseudo {pseudo} -sen {sen} -split {split} -splitExt {splitExt}".format(
            cmd = cmd, 
            pdir = cliParser.predir, 
            cut = cliParser.cut, 
            mcut = cliParser.mcut, 
            eps = cliParser.eps, 
            minPts = cliParser.minPts,
            pcut = cliParser.pcut,
            cpu = cliParser.cpu, 
            output = cliParser.fnOut, 
            bgd = cliParser.bgd,
            bgm = cliParser.bgm,
            pseudo = cliParser.pseudo, 
            sen = cliParser.sen,
            split = cliParser.split,
            splitExt = cliParser.splitExt,
        )
        logger.info(report)

        #check input
        if cliParser.predir == "":
            r = "No -d assigned! Return."
            logger.error(r)
            return
        if not os.path.exists(cliParser.predir):
            r = "%s not exists! Return." % cliParser.predir
            logger.error(r)
            return

        #check eps and minPts
        cliParser.eps = parseEps(cliParser.eps)
        if len(cliParser.eps) == 1 and cliParser.eps[0] == 0:
            logger.error("Input eps is 0. Return.")
            return
        cliParser.minPts = parseMinpts(cliParser.minPts)
        if len(cliParser.minPts) == 1 and cliParser.minPts[0] == 0:
            logger.error("Input minPts is 0. Return.")
            return

        if cliParser.bgd == "":
            callPeaks(
                cliParser.predir + "/petMeta.json",
                cliParser.fnOut,
                logger,
                eps=cliParser.eps,
                minPts=cliParser.minPts,
                pcut=cliParser.pcut,
                cpu=cliParser.cpu,
                pseudo=cliParser.pseudo,
                sen=cliParser.sen,
                cut=cliParser.cut,
                mcut=cliParser.mcut,
                split=cliParser.split,
                splitExt=cliParser.splitExt,
            )
        else:
            callPeaks(
                cliParser.predir + "/petMeta.json",
                cliParser.fnOut,
                logger,
                eps=cliParser.eps,
                minPts=cliParser.minPts,
                pcut=cliParser.pcut,
                cpu=cliParser.cpu,
                metabgf=cliParser.bgd + "/petMeta.json",
                bgm=cliParser.bgm,
                pseudo=cliParser.pseudo,
                sen=cliParser.sen,
                cut=cliParser.cut,
                mcut=cliParser.mcut,
                split=cliParser.split,
                splitExt=cliParser.splitExt,
            )

        end = datetime.now()
        logger.info("cLoops2 %s finished for %s. Used time: %s." %
                    (cmd, cliParser.fnOut, end - start) + "\n" * 3)

    #14. call loops
    if cmd == "callLoops":
        start = datetime.now()

        report = "Command: cLoops2 {} -d {} -eps {} -minPts {} -p {} -o {} -cut {} -mcut {} -filter {} -i {} -j {} -w {} -hic {} -max_cut {} -trans {} -emPair {}".format(
            cmd, 
            cliParser.predir, 
            cliParser.eps, 
            cliParser.minPts,
            cliParser.cpu, 
            cliParser.fnOut, 
            cliParser.cut,
            cliParser.mcut,
            cliParser.filterPETs, 
            cliParser.ucsc, 
            cliParser.juicebox,
            cliParser.washU, 
            cliParser.hic, 
            cliParser.max_cut,
            cliParser.trans,
            cliParser.emPair,
        )
        logger.info(report)

        #check input
        if cliParser.predir == "":
            r = "No -d assigned! Return."
            logger.error(r)
            return
        if not os.path.exists(cliParser.predir):
            r = "%s not exists! Return." % cliParser.predir
            logger.error(r)
            return

        #check eps and minPts
        cliParser.eps = parseEps(cliParser.eps,cliParser.emPair)
        if len(cliParser.eps) == 1 and cliParser.eps[0] == 0:
            logger.error("Input eps is 0. Return.")
            return
        cliParser.minPts = parseMinpts(cliParser.minPts,cliParser.emPair)
        if len(cliParser.minPts) == 1 and cliParser.minPts[0] == 0:
            logger.error("Input minPts is 0. Return.")
            return

        #check if loops file exits
        fout = cliParser.fnOut + "_loop.txt"
        if os.path.isfile(fout):
            r = "Output file %s exists, return." % fout
            logger.error(r)
            return

        #call cis loops
        callCisLoops(
            cliParser.predir,
            cliParser.fnOut,
            logger,
            eps=cliParser.eps,
            minPts=cliParser.minPts,
            cpu=cliParser.cpu,
            cut=cliParser.cut,
            mcut=cliParser.mcut,
            plot=cliParser.plot,
            max_cut=cliParser.max_cut,
            hic=cliParser.hic,
            filter=cliParser.filterPETs,
            ucsc=cliParser.ucsc,
            juicebox=cliParser.juicebox,
            washU=cliParser.washU,
            emPair=cliParser.emPair,
        )
        
        #call trans loops
        #check if loops file exits
        if cliParser.trans:
            fout = cliParser.fnOut + "_trans_loop.txt"
            if os.path.isfile(fout):
                r = "Output file %s exists, return." % fout
                logger.error(r)
                return

            #call loops
            callTransLoops(
                cliParser.predir,
                cliParser.fnOut,
                logger,
                eps=cliParser.eps,
                minPts=cliParser.minPts,
                cpu=cliParser.cpu,
                filter=cliParser.filterPETs,
                washU=cliParser.washU,
                juicebox=cliParser.juicebox,
            )


        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #15. call differentially enriched loops
    if cmd == "callDiffLoops":
        start = datetime.now()

        report = "Command: cLoops2 {} -tloop {} -cloop {} -td {} -cd {} -pcut {} -igp {} -noPCorr {} -fdr {} -o {} -p {} -j {} -w {} -customize {} -cacut {} -cmcut {} -vmin {} -vmax {} -cmap {}".format(
               cmd, 
               cliParser.tloop, 
               cliParser.cloop, 
               cliParser.tpred,
               cliParser.cpred, 
               cliParser.pcut, 
               cliParser.igp,
               cliParser.noPCorr,
               cliParser.fdr, 
               cliParser.fnOut,
               cliParser.cpu,
               cliParser.juicebox,
               cliParser.washU,
               cliParser.customize,
               cliParser.cacut,
               cliParser.cmcut,
               cliParser.vmin,
               cliParser.vmax,
               cliParser.cmap,
        )   
        logger.info(report)

        #check the file and directory
        if cliParser.tloop == "":
            r = "No -tloop assigned! Return."
            logger.error(r)
            return
        if not os.path.isfile(cliParser.tloop):
            r = "-tloop %s not exists! Return." % cliParser.tloop
            logger.error(r)
            return
        if cliParser.cloop == "":
            r = "No -cloop assigned! Return."
            logger.error(r)
            return
        if not os.path.isfile(cliParser.cloop):
            r = "-cloop %s not exists! Return." % cliParser.cloop
            logger.error(r)
            return
        if cliParser.tpred == "":
            r = "No -td assigned! Return."
            logger.error(r)
            return
        if not os.path.exists(cliParser.tpred):
            r = "-td %s not exists! Return." % cliParser.tpred
            logger.error(r)
            return
        if cliParser.cpred == "":
            r = "No -cd assigned! Return."
            logger.error(r)
            return
        if not os.path.exists(cliParser.cpred):
            r = "-cd %s not exists! Return." % cliParser.cpred
            logger.error(r)
            return
        if cliParser.customize:
            if cliParser.cacut == 0.0 and cliParser.cmcut == 0.0:
                r = "-customize option used, but -cacut and -cmcut not assigned! Return."
                logger.error(r)

        #run the call differentially enriched loops function
        callDiffLoops(
            cliParser.tloop,
            cliParser.cloop,
            cliParser.tpred,
            cliParser.cpred,
            cliParser.fnOut,
            cut=cliParser.cut,
            mcut=cliParser.mcut,
            cpu=cliParser.cpu,
            pcut=cliParser.pcut,
            igp=cliParser.igp,
            noPCorr=cliParser.noPCorr,
            fdrcut=cliParser.fdr,
            juicebox=cliParser.juicebox,
            washU=cliParser.washU,
            customize=cliParser.customize,
            cacut=cliParser.cacut,
            cmcut=cliParser.cmcut,
            vmin=cliParser.vmin,
            vmax=cliParser.vmax,
            cmap=cliParser.cmap,
        )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)
    
    #16. call domains
    if cmd == "callDomains":
        start = datetime.now()
        report = "Command: cLoops2 {} -d {} -cut {} -mcut {} -p {} -o {} -bs {} -ws {} -hic {}".format(
            cmd, 
            cliParser.predir, 
            cliParser.cut, 
            cliParser.mcut,
            cliParser.cpu, 
            cliParser.fnOut,
            cliParser.binSize, 
            cliParser.winSize, 
            cliParser.hic
        )
        logger.info(report)
        #check the input directory
        if cliParser.predir == "":
            logger.error("-d is required, return.")
            return
        cliParser.binSize = parseEps(cliParser.binSize)
        cliParser.winSize = parseEps(cliParser.winSize)
        if len(cliParser.binSize) == 1 and cliParser.binSize[0] == 0:
            logger.error(
                "Input -bs is 0, please input with possible reasonable resoultion. Return."
            )
            return

        #do the job
        callDomains(
            cliParser.predir + "/petMeta.json",
            cliParser.fnOut,
            logger,
            bs=cliParser.binSize,
            ws=cliParser.winSize,
            cut=cliParser.cut,
            mcut=cliParser.mcut,
            cpu=cliParser.cpu,
            hic=cliParser.hic
            )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #17. plot matrix, peaks, loops, domains
    if cmd == "plot":
        start = datetime.now()

        report = "Command cLoops2 {cmd} -f {f} -o {output} -chrom {chrom} -start {start} -end {end} -bs {r} -cut {cut} -mcut {mcut} -log {log} -m {method} -corr {corr} -triu {triu} -norm {norm} -bws {bws} -bwvs {bwvs} -bwcs {bwcs} -beds {beds} -gtf {gtf} -1D {oneD} -1Dv {oneDv} -loops {floop} -domains {fdomain} -vmin {vmin} -vmax {vmax} -virtual4C {virtual4C} -view_start {viewStart} -view_end {viewEnd} -arch {arch} -aw {aw} -ac {ac} -aa {aa} -scatter {scatter} -ss {ss} -sa {sa} -sc {sc} -eig {eig} -eig_r {eig_r} -figWidth {figWidth}".format(
            cmd=cmd,
            f=cliParser.fixy,
            output=cliParser.fnOut,
            chrom=cliParser.chrom,
            start=cliParser.start,
            end=cliParser.end,
            r=cliParser.binSize,
            cut=cliParser.cut,
            mcut=cliParser.mcut,
            log=cliParser.log,
            method=cliParser.method,
            corr=cliParser.corr,
            triu=cliParser.triu,
            norm=cliParser.norm,
            bws=cliParser.bws,
            bwvs=cliParser.bwvs,
            bwcs=cliParser.bwcs,
            beds=cliParser.beds,
            gtf=cliParser.gtf,
            oneD=cliParser.oneD,
            oneDv=cliParser.oneDv,
            floop=cliParser.floop,
            fdomain=cliParser.fdomain,
            vmin=cliParser.vmin,
            vmax=cliParser.vmax,
            virtual4C=cliParser.virtual4C,
            viewStart=cliParser.viewStart,
            viewEnd=cliParser.viewEnd,
            arch=cliParser.arch,
            aw=cliParser.aw,
            ac=cliParser.ac,
            aa=cliParser.aa,
            scatter=cliParser.scatter,
            ss=cliParser.ss,
            sa=cliParser.sa,
            sc=cliParser.sc,
            eig=cliParser.eig,
            eig_r=cliParser.eig_r,
            figWidth=cliParser.figWidth,
        )
        logger.info(report)

        #check everything
        if cliParser.binSize <=0:
            logger.error("ERROR! -bs %s is <=0!"%cliParser.binSize)
            return
        if cliParser.start != -1 and cliParser.end != -1 and cliParser.end < cliParser.start:
            logger.error("ERROR! end %s is smaller than %s start." %
                         (cliParser.end, cliParser.start))
            return
        if cliParser.bws != "":
            bws = cliParser.bws.split(",")
            bws = [bw for bw in bws if os.path.isfile(bw)]
        else:
            bws = []
        if cliParser.bwvs !="":
            bwvs = cliParser.bwvs.split(";")
            if len(bwvs)!=len(bws):
                logger.error("ERROR! Length of input -bws %s and -bwvs not equal" %(cliParser.bws, cliParser.bwvs))
                return
        if cliParser.bwcs != "":
            bwcs = cliParser.bwcs.split(",")
            if len(bwcs)!=len(bws):
                logger.error("ERROR! Length of input -bws %s and -bwcs not equal" %(cliParser.bws, cliParser.bwcs))
                return
        if cliParser.beds != "":
            beds = cliParser.beds.split(",")
            beds = [bed for bed in beds if os.path.isfile(bed)]
        else:
            beds = []
        if cliParser.floop != "":
            try:
                loops = parseTxt2Loops(cliParser.floop)
            except:
                logger.error(
                    "Error in parsing file %s as loops, please check." %
                    (cliParser.floop))
                return
        else:
            loops = None
        if cliParser.fdomain != "":
            if not os.path.isfile(cliParser.fdomain):
                logger.error("Domain file %s not exists, please check." %
                             (cliParser.fdomain))
                return
        if cliParser.virtual4C:
            if cliParser.end !=-1 and not cliParser.start <= cliParser.start <=cliParser.viewStart \
                <= cliParser.viewEnd <= cliParser.end:
                logger.error("Please check -view_start and -view_end, should be start<=view_start<=view_end<=end")
                return 
        if cliParser.gtf !="":
            if not os.path.isfile(cliParser.gtf):
                logger.error("Gene track of -gtf %s not exists, please check."%(cliParser.gtf))
                return

        if cliParser.fixy != "":
            if not os.path.isfile(cliParser.fixy):
                logger.info("ERROR! -f assigned %s but not exists. Return."%cliParser.fixy)
                return
            #plot heatmap
            if cliParser.arch == False and cliParser.scatter == False:
                plotMatHeatmap(
                    cliParser.fixy,
                    cliParser.fnOut,
                    start=cliParser.start,
                    end=cliParser.end,
                    res=cliParser.binSize,
                    cut=cliParser.cut,
                    mcut=cliParser.mcut,
                    log=cliParser.log,
                    method=cliParser.method,
                    oneD=cliParser.oneD,
                    oneDv=cliParser.oneDv,
                    corr=cliParser.corr,
                    triu=cliParser.triu,
                    norm=cliParser.norm,
                    bws=bws,
                    bwvs=cliParser.bwvs,
                    bwcs=cliParser.bwcs,
                    beds=beds,
                    gtf=cliParser.gtf,
                    loops=loops,
                    domains=cliParser.fdomain,
                    virtual4C=cliParser.virtual4C,
                    viewStart=cliParser.viewStart,
                    viewEnd=cliParser.viewEnd,
                    vmin=cliParser.vmin,
                    vmax=cliParser.vmax,
                    eig=cliParser.eig,
                    eig_r=cliParser.eig_r,
                    width=cliParser.figWidth,
                )
            else:
                if cliParser.arch:
                    plotPETsArches( 
                        cliParser.fixy,
                        cliParser.fnOut, 
                        start=cliParser.start, 
                        end=cliParser.end, 
                        cut=cliParser.cut, 
                        mcut=cliParser.mcut, 
                        oneD=cliParser.oneD, 
                        oneDv=cliParser.oneDv,
                        bws=bws, 
                        bwvs=cliParser.bwvs,
                        bwcs=cliParser.bwcs,
                        beds=beds, 
                        loops=loops, 
                        gtf=cliParser.gtf, 
                        aw=cliParser.aw,
                        ac=cliParser.ac,
                        aa=cliParser.aa,
                        width=cliParser.figWidth,
                    )
                if cliParser.scatter:
                    plotPETsScatter( 
                        cliParser.fixy,
                        cliParser.fnOut, 
                        start=cliParser.start, 
                        end=cliParser.end, 
                        cut=cliParser.cut, 
                        mcut=cliParser.mcut, 
                        oneD=cliParser.oneD, 
                        oneDv=cliParser.oneDv,
                        bws=bws, 
                        bwvs=cliParser.bwvs,
                        bwcs=cliParser.bwcs,
                        beds=beds, 
                        loops=loops, 
                        gtf=cliParser.gtf, 
                        ss=cliParser.ss,
                        sc=cliParser.sc,
                        sa=cliParser.sa,
                        triu=cliParser.triu,
                        width=cliParser.figWidth,
                    )
        else:
            if cliParser.chrom !="":
                plotProfiles( 
                    cliParser.fnOut, 
                    chrom=cliParser.chrom,
                    start=cliParser.start, 
                    end=cliParser.end, 
                    bws=bws, 
                    bwvs=cliParser.bwvs,
                    bwcs=cliParser.bwcs,
                    beds=beds, 
                    loops=loops, 
                    gtf=cliParser.gtf, 
                    width=cliParser.figWidth,
                )
            else:
                logger.error("ERROR! -f %s not exists and -chrom %s not assigned"%(cliParser.fixy,cliParser.chrom))
                return 
 
        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)
    
    #18. montage analysis 
    if cmd == "montage":
        start = datetime.now()

        report = "Command cLoops2 {cmd} -f {f} -bed {bed} -ext {ext} -o {output} -cut {cut} -mcut {mcut} -simple {simple} -viewPoint {viewPoint} -vmin {vmin} -vmax {vmax} -ppmw {ppmw} -aw {aw} -no1D {noOneD}".format(
            cmd=cmd,
            f=cliParser.fixy,
            bed=cliParser.bed,
            output=cliParser.fnOut,
            cut=cliParser.cut,
            mcut=cliParser.mcut,
            ext=cliParser.ext,
            simple=cliParser.simple,
            viewPoint=cliParser.viewPoint, 
            vmin=cliParser.vmin,
            vmax=cliParser.vmax,
            ppmw=cliParser.ppmw,
            aw=cliParser.aw,
            noOneD=cliParser.noOneD,
        )
        logger.info(report)
        
        if cliParser.noOneD:
            oneD = False
        else:
            oneD = True

        #run analysis
        montage(
            cliParser.fixy,
            cliParser.bed,
            cliParser.fnOut,
            ext=cliParser.ext,
            cut=cliParser.cut,
            mcut=cliParser.mcut,
            simple=cliParser.simple,
            viewPoint=cliParser.viewPoint,
            vmin=cliParser.vmin,
            vmax=cliParser.vmax,
            ppmw=cliParser.ppmw,
            aw=cliParser.aw,
            oneD=oneD,
        )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #19. aggregated analysis and plot
    if cmd == "agg":
        start = datetime.now()

        report = "Command cLoops2 {cmd} -d {predir} -o {output} -cut {cut} -mcut {mcut} -p {cpu} -skipZeros {skipZeros} -peaks {peakf} -peak_ext {peak_ext} -peak_bins {peak_bins} -peak_norm {peak_norm} -loops {loopf} -loop_ext {loop_ext} -loop_cut {loop_cut} -loop_norm {loop_norm} -loop_vmin {loop_vmin} -loop_vmax {loop_vmax} -viewPoints {viewPointF} -viewPointUp {viewPointUp} -viewPointDown {viewPointDown} -viewPointBs {viewPointBs} -viewPoint_norm {viewPoint_norm} -viewPoint_vmin {viewPoint_vmin} -viewPoint_vmax {viewPoint_vmax} -twoAnchors {twoAnchorF} -twoAnchor_ext {twoAnchor_ext} -twoAnchor_vmin {twoAnchor_vmin} -twoAnchor_vmax {twoAnchor_vmax} -domains {domainf} -domain_ext {domain_ext} -domain_vmin {dvmin} -domain_vmax {dvmax} -bws {bws} -1D {oneD}".format(
            cmd=cmd,
            predir=cliParser.predir,
            output=cliParser.fnOut,
            cut=cliParser.cut,
            mcut=cliParser.mcut,
            cpu=cliParser.cpu,
            peakf=cliParser.peakf,
            peak_ext=cliParser.peak_ext,
            peak_bins=cliParser.peak_bins,
            peak_norm=cliParser.peak_norm,
            loopf=cliParser.loopf,
            loop_ext=cliParser.loop_ext,
            loop_cut=cliParser.loop_cut,
            loop_norm=cliParser.loop_norm,
            loop_vmin=cliParser.loop_vmin,
            loop_vmax=cliParser.loop_vmax,
            viewPointF=cliParser.viewPointF,
            viewPointUp=cliParser.viewPointUp,
            viewPointDown=cliParser.viewPointDown,
            viewPointBs=cliParser.viewPointBs,
            viewPoint_norm=cliParser.viewPoint_norm,
            viewPoint_vmin=cliParser.viewPoint_vmin,
            viewPoint_vmax=cliParser.viewPoint_vmax,
            twoAnchorF=cliParser.twoAnchorsF,
            twoAnchor_ext=cliParser.twoAnchor_ext,
            twoAnchor_vmin=cliParser.twoAnchor_vmin,
            twoAnchor_vmax=cliParser.twoAnchor_vmax,
            domainf=cliParser.domainf,
            domain_ext=cliParser.domain_ext,
            dvmin=cliParser.domain_vmin,
            dvmax=cliParser.domain_vmax,
            bws=cliParser.bws,
            oneD=cliParser.oneD,
            skipZeros=cliParser.skipZeros,
        )
        logger.info(report)

        #check input data directory
        if cliParser.predir == "":
            r = "No -d assigned! Return."
            logger.error(r)
            return
        if not os.path.isdir(cliParser.predir):
            r = "%s not exists! Return."%(cliParser.predir)
            logger.error(r)
            return
        if cliParser.peakf == "" and cliParser.loopf == "" and cliParser.domainf == "" and cliParser.viewPointF == "" and cliParser.twoAnchorsF == "":
            logger.error("No -peaks, -loops, -viewPoints, -domains or -twoAnchors assigned, nothing to analyze, return.")
            return

        #run aggregation peaks analysis
        if cliParser.peakf != "":
            if os.path.isfile(cliParser.peakf):
                aggPeaks(
                    cliParser.predir,
                    cliParser.peakf,
                    cliParser.fnOut,
                    logger,
                    ext=cliParser.peak_ext,
                    bins=cliParser.peak_bins,
                    cut=cliParser.cut,
                    mcut=cliParser.mcut,
                    cpu=cliParser.cpu,
                    skipZeros=cliParser.skipZeros,
                    norm=cliParser.peak_norm,
                )
            else:
                logger.error("%s not exists! Return."%cliParser.peakf)
        
        #run aggregation loops analysis
        if cliParser.loopf != "":
            if os.path.isfile(cliParser.loopf):
                aggLoops(
                    cliParser.predir,
                    cliParser.loopf,
                    cliParser.fnOut,
                    logger,
                    bws=cliParser.bws,
                    ext=cliParser.loop_ext,
                    cut=cliParser.cut,
                    mcut=cliParser.mcut,
                    lcut=cliParser.loop_cut,
                    cpu=cliParser.cpu,
                    skipZeros=cliParser.skipZeros,
                    norm=cliParser.loop_norm,
                    oneD=cliParser.oneD,
                    vmin=cliParser.loop_vmin,
                    vmax=cliParser.loop_vmax,
                )
            else:
                logger.error("%s not exists! Return."%cliParser.loopf)

        #run aggregation virtual 4C view points analysis
        if cliParser.viewPointF != "":
            if os.path.isfile(cliParser.viewPointF):
                aggViewPoints(
                    cliParser.predir,
                    cliParser.viewPointF,
                    cliParser.fnOut,
                    logger,
                    bws=cliParser.bws,
                    upExt=cliParser.viewPointUp,
                    downExt=cliParser.viewPointDown,
                    bs=cliParser.viewPointBs,
                    cut=cliParser.cut,
                    mcut=cliParser.mcut,
                    cpu=cliParser.cpu,
                    skipZeros=cliParser.skipZeros,
                    oneD=cliParser.oneD,
                    norm=cliParser.viewPoint_norm,
                    vmin=cliParser.viewPoint_vmin,
                    vmax=cliParser.viewPoint_vmax,
                )
            else:
                logger.error("%s not exists! Return."%cliParser.viewPointF)

        #run aggregation domains analysis
        if cliParser.domainf != "":
            if os.path.isfile(cliParser.domainf):
                aggDomains(
                    cliParser.predir,
                    cliParser.domainf,
                    cliParser.fnOut,
                    logger,
                    bws=cliParser.bws,
                    ext=cliParser.domain_ext,
                    cut=cliParser.cut,
                    mcut=cliParser.mcut,
                    cpu=cliParser.cpu,
                    vmin=cliParser.domain_vmin,
                    vmax=cliParser.domain_vmax,
                    skipZeros=cliParser.skipZeros,
                    oneD=cliParser.oneD,
                )
            else:
                logger.error("%s not exists! Return."%cliParser.domainf)

        #run aggregation domains analysis
        if cliParser.twoAnchorsF != "":
            if os.path.isfile(cliParser.twoAnchorsF):
                aggTwoAnchors(
                    cliParser.predir,
                    cliParser.twoAnchorsF,
                    cliParser.fnOut,
                    logger,
                    ext=cliParser.twoAnchor_ext,
                    bws=cliParser.bws,
                    cut=cliParser.cut,
                    mcut=cliParser.mcut,
                    cpu=cliParser.cpu,
                    skipZeros=cliParser.skipZeros,
                    oneD=cliParser.oneD,
                    vmin=cliParser.twoAnchor_vmin,
                    vmax=cliParser.twoAnchor_vmax,
                )
            else:
                logger.error("%s not exists! Return."%cliParser.twoAnchorsF)

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)
    
    #20. quantify peaks,loops and domains
    if cmd == "quant":
        start = datetime.now()

        report = "Command cLoops2 {cmd} -d {predir} -o {output} -cut {cut} -mcut {mcut} -p {cpu} -peaks {peakf} -loops {loopf} -domains {domainf} -domain_bs {bs} -domain_ws {ws} -domain_bdg {dbdg}".format(
            cmd=cmd,
            predir=cliParser.predir,
            output=cliParser.fnOut,
            cut=cliParser.cut,
            mcut=cliParser.mcut,
            cpu=cliParser.cpu,
            peakf=cliParser.peakf,
            loopf=cliParser.loopf,
            domainf=cliParser.domainf,
            bs=cliParser.domainBinSize,
            ws=cliParser.domainWinSize,
            dbdg=cliParser.domainBdg,
        )
        logger.info(report)

        #check input data directory
        if cliParser.predir == "":
            r = "No -d assigned! Return."
            logger.error(r)
            return
        if not os.path.isdir(cliParser.predir):
            r = "%s not exists! Return."%(cliParser.predir)
            logger.error(r)
            return
        if cliParser.peakf == "" and cliParser.loopf == "" and cliParser.domainf == "":
            logger.error("No -peaks, -loops, -domains assigned, nothing to analyze, return.")
            return

        #run aggregation peaks analysis
        if cliParser.peakf != "":
            if os.path.isfile(cliParser.peakf):
                quantPeaks(
                    cliParser.predir,
                    cliParser.peakf,
                    cliParser.fnOut,
                    logger,
                    cut=cliParser.cut,
                    mcut=cliParser.mcut,
                    cpu=cliParser.cpu,
                )
            else:
                logger.error("%s not exists! Return."%cliParser.peakf)
        if cliParser.loopf != "":
            if os.path.isfile(cliParser.loopf):
                quantLoops(
                    cliParser.predir,
                    cliParser.loopf,
                    cliParser.fnOut,
                    logger,
                    cut=cliParser.cut,
                    mcut=cliParser.mcut,
                    cpu=cliParser.cpu,
                )
            else:
                logger.error("%s not exists! Return."%cliParser.loopf)
        if cliParser.domainf != "":
            if os.path.isfile(cliParser.domainf):
                quantDomains(
                    cliParser.predir,
                    cliParser.domainf,
                    cliParser.fnOut,
                    logger,
                    bs=cliParser.domainBinSize,
                    ws=cliParser.domainWinSize,
                    cut=cliParser.cut,
                    mcut=cliParser.mcut,
                    cpu=cliParser.cpu,
                    bdg=cliParser.domainBdg
                )
            else:
                logger.error("%s not exists! Return."%cliParser.domainf)

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #21. analyze loops 
    if cmd == "anaLoops":
        start = datetime.now()

        report = "Command cLoops2 {cmd} -loops {floop} -o {output} -gtf {gtf} -tid {tid} -p {cpu} -pdis {pdis} -net {net} -gap {gap}".format(
            cmd=cmd,
            floop=cliParser.floop,
            output=cliParser.fnOut,
            gtf=cliParser.gtf,
            tid=cliParser.tid,
            cpu=cliParser.cpu,
            pdis=cliParser.pdis,
            net=cliParser.net,
            gap=cliParser.gap,
        )
        logger.info(report)

        if cliParser.floop == "" or not os.path.isfile( cliParser.floop ):
            logger.error("%s not exists! Return."%cliParser.floop)
            return
        
        #do the job
        anaLoops(
            cliParser.floop,
            cliParser.fnOut,
            gtf=cliParser.gtf,
            tid=cliParser.tid,
            cpu=cliParser.cpu,
            pdis=cliParser.pdis,
            net=cliParser.net,
            gap=cliParser.gap,
        )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)

    #22. find loops/regions target genes 
    if cmd == "findTargets":
        start = datetime.now()

        report = "Command cLoops2 {cmd} -net {net} -tg {tg} -bed {bed} -o {output} ".format(
            cmd=cmd,
            net=cliParser.fnet,
            tg=cliParser.ftg,
            bed=cliParser.fbed,
            output=cliParser.fnOut,
        )
        logger.info(report)
        
        if cliParser.fbed !="" and not os.path.isfile(cliParser.fbed):
            logger.error("%s not exists! Return."%cliParser.fbed)
            return

        #do the job
        findTargets(
            cliParser.fnet,
            cliParser.ftg,
            cliParser.fnOut,
            fbed=cliParser.fbed,
        )

        end = datetime.now()
        logger.info("cLoops2 %s finished. Used time: %s." %
                    (cmd, end - start) + "\n" * 3)



    
if __name__ == "__main__":
    main()
