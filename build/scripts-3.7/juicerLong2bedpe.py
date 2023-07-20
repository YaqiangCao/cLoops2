#!/usr/bin/env pyhthon
#--coding:utf-8--
"""
"""

import argparse, gzip, os, sys
from datetime import datetime
from argparse import RawTextHelpFormatter



def help():
    description = """
        Convert Juicer long format file to to BEDPE file as input of cLoops2. 
        Example:
        juicerLong2bedpe.py -f test.allValidPairs -o test
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-f',
                        dest="fin",
                        type=str,
                        help="Input file name, required.")
    parser.add_argument('-o',
                        dest="fout",
                        required=True,
                        type=str,
                        help="Output file name, required.")
    op = parser.parse_args()
    return op



def long2bedpe(fin, fout, ext=50):
    with open(fout, "w") as fo:
        for line in open(fin):
            if line.startswith("#"):
                continue
            line = line.split("\n")[0].split()
            nline = [
                line[1],
                max(0,
                    int(line[2]) - ext),
                int(line[2]) + ext,  #pet 1 
                line[5],
                max(0,
                    int(line[6]) - ext),
                int(line[6]) + ext,  #pet 2
                ".",
                ".",
                "+",
                "+"  #other infor
            ]
            if line[0] != "0":
                nline[-2] = "-"
            if line[4] != "0":
                nline[-1] = "-"
            fo.write("\t".join(list(map(str, nline))) + "\n")



def main():
    op = help()
    if not os.path.isfile(op.fin):
        sys.stderr.write("Error: input file %s not exists!\n" % op.fin)
    if os.path.isfile(op.fout):
        sys.stderr.write("Error: output file %s exists! \n" % op.fout)
    long2bedpe(op.fin, op.fout)


if __name__ == '__main__':
    start_time = datetime.now()
    main()
    usedtime = datetime.now() - start_time
    sys.stderr.write("Process finished. Used CPU time: %s Bye!\n" % usedtime)
