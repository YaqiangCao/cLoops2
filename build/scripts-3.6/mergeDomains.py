#!python
#--coding:utf-8 --
"""
Merge domains from multiple results of different resolutions. 
"""

#sys
import os
import sys
import argparse
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library

#cLoops2
from cLoops2.ds import Domain
from cLoops2.io import doms2bed



def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Merge domains from multiple resolutions.

        Example:
        mergeDomains.py -fs 5k_tad.bed,10k_tad.bed,25k_tad.bed -o all -r 0.9
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-fs",
        dest="fins",
        required=True,
        type=str,
        help=
        "The domains in BED format from differnet resolutions to be merged.\n"\
        "BED files should be input with decreased resolutions, for example\n"\
        "5k.bed,10k.bed,25k.bed. Multiple files seperated by comma."
    )
    parser.add_argument(
        '-r',
        dest="lrcut",
        required=False,
        type=float,
        default=0.9,
        help="Ratio of overlapped domains to be merged. Default is 0.9."
    )
    parser.add_argument(
        '-o',
        dest="fout",
        required=True,
        type=str,
        help="Output file name, required."
    )
    op = parser.parse_args()
    return op



def compDoms(doma,domb,lrcut=0.9):
    """
    Compare if is quite close same domains.
    If quite close, whether use doma to replace domb.
    """
    if doma.chrom != domb.chrom:
        return False
    #overlapped domains
    if domb.start <= doma.start <= domb.end or domb.start <= doma.end <= domb.end or doma.start <= domb.start <= doma.end or doma.start <= domb.end <= doma.end:
        start = max(doma.start,domb.start)
        end = min(doma.end,domb.end)
        length = max(doma.length,domb.length)
        if (end-start)/length > lrcut:
            return True
    return False
 

def combineDoms(doms, doms2, lrcut=0.9):
    """
    Combine domains.
    """
    #doms binsize is bigger than doms2
    for key in doms2.keys():
        if key not in doms:
            doms[key] = doms2[key]
        else:
            #add no overlapped
            for doma in doms2[key]:
                flag = False
                for i, domb in enumerate(doms[key]):
                    flag2 = compDoms(doma,domb,lrcut) 
                    #if highly overlapped and similar, use the higher resolution/smaller domain to replace biggger one
                    if flag2:
                        flag = True
                        doms[key][i] = doma #replace 
                        break
                    else:
                        continue
                #no overlapped or almost the same
                if flag == False:
                    doms[key].append(doma)
    return doms


def readDoms(f):
    doms = {}
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        if len(line) == 3:
            did = "|".join(line)
        if len(line) >=4:
            if line[4].strip() == "":
                did = "|".join(line)
            else:
                did = line[4]
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        dom = Domain()
        dom.id = did
        dom.chrom = chrom
        dom.start = start
        dom.end = end
        dom.length = end - start
        if chrom not in doms:
            doms[chrom] = []
        doms[chrom].append( dom )
    return doms
        
                

def main():
    op = help()
    
    fs = op.fins.split(",")
    for f in fs:
        if not os.path.isfile(f):
            print("%s not exists, return!"%f)
    
    doms = readDoms( fs[0] )
    for f in fs[1:]:
        doms = combineDoms( doms, readDoms(f), lrcut=op.lrcut)
    
    with open(op.fout+"_mergedDomains.bed","w") as fo:
        for c,v in doms.items():
            for dom in v:
                line = list( map(str, [dom.chrom, dom.start, dom.end, dom.id]) )
                fo.write( "\t".join(line) + "\n")
    with open(op.fout+"_mergedDomains_juicebox.bedpe","w") as fo:
        head = "chr1\tx1\tx2\tchr2\ty1\ty2\tname\n"
        fo.write(head)
        for c,v in doms.items():
            for dom in v:
                line = list( map(str, [dom.chrom, dom.start, dom.end,dom.chrom,dom.start,dom.end, dom.id]) )
                fo.write( "\t".join(line) + "\n")
    

if __name__ == "__main__":
    main()
