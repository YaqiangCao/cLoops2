#!/usr/bin/env python
#--coding:utf-8--
"""
hicpro2bedpe.py
Convert HiC-Pro output validate pairs to cLoops input bedpe file.
"""
#sys
import os
import re
import sys
import gzip
import argparse
from glob import glob
from datetime import datetime
from argparse import RawTextHelpFormatter

#3rd

#cLoops2
from cLoops2.utils import cFlush


def help():
    description = """
        Convert HiC-Pro allValidParis to BEDPE file as input of cLoops2. 
        Example:
        hicpro2bedpe.py -f test.allValidPairs -o test
        """

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        '-f',
        dest="fin",
        required=True,
        type=str,
        help="HiC-Pro allValidPairs file, can be .allValidParis.gz. ")
    parser.add_argument('-o',
                        dest="out",
                        required=True,
                        type=str,
                        help="Output file prefix.")
    parser.add_argument('-ext',
                        dest="ext",
                        required=False,
                        type=int,
                        default=50,
                        help="Extension from center of read, default is 50. ")
    op = parser.parse_args()
    return op


def pairs2bedpe(f_hicpro, f_out,ext=50):
    """
    Converting HiC-Pro output allValidPairs to bedpe file.
    """
    with gzip.open(f_out, 'wt') as f_bedpe:
        if f_hicpro.endswith('.gz'):
            #f_pair = gzip.open(f_hicpro) #python2
            f_pair = gzip.open(f_hicpro, 'rt')  #python3
        else:
            f_pair = open(f_hicpro)
        for i, line in enumerate(f_pair):
            if i % 100000 == 0:
                cFlush("%s PETs processed from %s" % (i, f_hicpro))
            line = line.strip().split('\t')
            #if the position is middle of reads
            #petA = [line[1], int(line[2])-ext, int(line[2])+ext]
            #petB = [line[4], int(line[5])-ext, int(line[5])+ext]
            #if the position is 5 end of reads
            if line[3] == "+":
                petA = [line[1], int(line[2]), int(line[2])+ext]
            else:
                petA = [line[1], int(line[2])-ext, int(line[2])]
            if line[6] == "+":
                petB = [line[4], int(line[5]), int(line[5])+ext]
            else:
                petB = [line[4], int(line[5])-ext, int(line[5])]

            newline = [
                petA[0], petA[1], petA[2], petB[0], petB[1], petB[2], line[0],
                '.', line[3], line[6]
            ]
            f_bedpe.write("\t".join(map(str, newline)) + "\n")
        f_pair.close()


def main():
    op = help()
    if not os.path.isfile(op.fin):
        sys.stderr.write("Error: input file %s not exists.\n" % op.fin)
        return
    bedpe_file = re.sub(r'_allValidPairs(.gz)?$', '', op.out)
    bedpe_file = bedpe_file + '.bedpe.gz'
    if os.path.isfile(op.out):
        sys.stderr.write("Error: output file %s exists.\n" % bedpe_file)
        return
    pairs2bedpe(op.fin, bedpe_file,ext=op.ext)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    usedtime = datetime.now() - start_time
    sys.stderr.write("Process finished. Used CPU time: %s Bye!\n" % usedtime)
