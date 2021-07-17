#!python
#--coding:utf-8--
"""
getBedpeFBed.py
Transfering single-end BED file to paired-end BEDPE file as input of cLoops2 .
"""
#systematic library
import os, time, gzip, argparse, sys
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd library

#cLoops2
from cLoops2.ds import PET
from cLoops2.utils import cFlush


def help():
    description = """
        Transfering single-end BED file to paired-end BEDPE file as input of
        cLoops2 for furthur analysis.
        The 6th column of the BED file of strand information is used to extend 
        the fragments.
        If no strand information available, default treat it as + strand
        Example:
        getBedpeFBed.py -f a.bed.gz -o a
        """

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        '-f',
        dest="fin",
        required=True,
        type=str,
        help=
        "Input bed files, or .bed.gz files. "
    )
    parser.add_argument('-o',
                        dest="out",
                        required=True,
                        type=str,
                        help="Output file prefix.")
    parser.add_argument(
        '-ext',
        dest="ext",
        required=False,
        default=150,
        type=int,
        help=
        "The expect fragment length of the bed file to extend from 5' to 3', default is 150."
    )

    op = parser.parse_args()
    return op


def bed2bedpe(fin, fout, ext=150):
    """
    Extend the BED file to BEDPE file according to expected fragment size.
    """
    if fin.endswith(".gz"):
        fino = gzip.open(fin, "rt")
    else:
        fino = open(fin)
    if fout.endswith(".gz"):
        fo = gzip.open(fout, "wt")
    else:
        fo = open(fout, "w")
    for i, line in enumerate(fino):
        if i % 10000 == 0:
            cFlush("%s read from %s" % (i,fin))
        line = line.split("\n")[0].split('\t')
        if len(line) < 6:  #no strand information
            nline = [
                line[0], line[1], line[2], line[0],
                int(line[1]) + ext,
                int(line[2]) + ext, ".", "44", "+", "-"
            ]
        elif line[5] == "+":
            nline = [
                line[0], line[1], line[2], line[0],
                int(line[1]) + ext,
                int(line[2]) + ext, ".", "44", "+", "-"
            ]
        else:
            nline = [
                line[0],
                max(0, int(line[1])),
                max(0,
                    int(line[2]) - ext), line[0], line[1], line[2], ".", "44",
                "+", "-"
            ]
        nline = "\t".join(list(map(str, nline))) + "\n"
        fo.write(nline)
    fino.close()
    fo.close()


def main():
    op = help()
    bed2bedpe(op.fin, op.out+".bedpe.gz", ext=op.ext)


if __name__ == "__main__":
    start_time = datetime.now()
    main()
    usedtime = datetime.now() - start_time
    sys.stderr.write("Process finished. Used CPU time: %s Bye!\n" % usedtime)
