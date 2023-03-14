#!/home/caoy7/anaconda2/envs/cLoops2/bin/python
#--coding:utf-8 --
"""
compareComp.py
cLoops2 compareComp.py compare compartment PC1 values based on Mahalanobis distance and annotate the changed bins.

"""

__date__ = "2023-03-09"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys library
import argparse
from argparse import RawTextHelpFormatter

#3rd library
import pylab
import HTSeq
import pandas as pd
import numpy as np

#cLoops2
from cLoops2.stat import twoPassesMDTest
from cLoops2.settings import *


def help():
    """
    Create the command line interface for the script.
    """
    description = """
        Pair-wisely compare compartments PC1 based on two-passes Mahalanobis distance. 

        Example:
        compareComp.py -a young_pc1.bdg -b old_pc1.bdg -o youngVsOld -na Young -b Old -pcut 0.01 -gtf mm10.gtf
        """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-a",
        dest="fa",
        required=True,
        type=str,
        help=
        "Input bedGraph file for the first sample containing the PC1 scores from comparment analysis. Should be aligned for A/B considering gene density or CpG islands contents. Better to be a control such as young or wild-type sample."
    )
    parser.add_argument(
        "-b",
        dest="fb",
        required=True,
        type=str,
        help=
        "Input bedGraph file for the second sample containing the PC1 scores from comparment analysis."
    )
    parser.add_argument(
        "-na",
        dest="na",
        type=str,
        required=True,
        help=
        "Name of first sample, will be shown in the output figure. Only use alphabet and numbers."
    )
    parser.add_argument(
        "-nb",
        dest="nb",
        type=str,
        required=True,
        help="Name of second sample, will be shown in the output figure.")
    parser.add_argument(
        "-gtf",
        dest="gtf",
        default="",
        required=True,
        type=str,
        help=
        "GTF file annotation for genes. Significant flip/switch overlapped genes will be reported based on the gene annotation file."
    )
    parser.add_argument(
        "-pcut",
        dest="pcut",
        type=float,
        default=0.01,
        help=
        "Chi-Square p-value cutoff for detecting siginficant different compartment, default is 0.01."
    )

    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
    op = parser.parse_args()
    return op


def stichBins2Comp(rs):
    """
    Stich bins to compartment
    """
    comps = []
    i = 0
    while i < len(rs):
        j = i + 1
        p = i
        while j < len(rs):
            if rs[j][-1] * rs[i][-1] > 0:
                p = j
                j += 1
            else:
                break
        if rs[i][-1] > 0:
            flag = "compartmentA"
        else:
            flag = "compartmentB"
        s = rs[i][1]
        e = rs[p][2]
        vs = [rs[t][-1] for t in range(i, j)]
        v = np.mean(vs)
        comps.append([
            rs[i][0], s, e,
            rs[i][0] + ":" + str(s) + "-" + str(e) + "|" + flag,
            str(v)
        ])
        i = j
    return comps


def getCompartment(bdgf, fout):
    """
    Get compartment according to PC1 values.
    """
    ds = {}
    for line in open(bdgf):
        line = line.split("\n")[0].split("\t")
        if line[0] not in ds:
            ds[line[0]] = []
        line[-1] = float(line[-1])
        ds[line[0]].append(line)
    #stich domains
    cs = list(ds.keys())
    cs.sort()
    with open(fout, "w") as fo:
        for c in cs:
            rs = stichBins2Comp(ds[c])
            for r in rs:
                fo.write("\t".join(r) + "\n")


def readBdg(f):
    """
    Read bedGraph file for PC1.
    """
    ds = {}  #pandas series
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        k = line[0] + ":" + line[1] + "-" + line[2]
        v = float(line[-1])
        ds[k] = v
    ds = pd.Series(ds)
    return ds


def getFlips(sa, sb, ps, inds):
    """
    Get compartment flips according to PC1 and p-values. 

    @param sa: pd.Series, PC1 values for sample A
    @param sb: pd.Series, PC1 values for sample B
    @param ps: pd.Series, p-values
    @param inds: pd.Series.idnex, p-value select index

    @return a2b: pd.Series.idnex, A to B flips 
    @return b2a: pd.Series.idnex, B to A flips
    @return a2a: pd.Series.idnex, A to A changes 
    @return b2b: pd.Series.index, B to B changes
    """
    #AtoB flips
    ta = sa[sa > 0].index
    tb = sb[sb < 0].index
    a2b = ta.intersection(tb).intersection(inds)
    a2b = ps[a2b].sort_values(inplace=False, ascending=True).index

    #BtoA flips
    ta = sa[sa < 0].index
    tb = sb[sb > 0].index
    b2a = ta.intersection(tb).intersection(inds)
    b2a = ps[b2a].sort_values(inplace=False, ascending=True).index

    #AtoA
    ta = sa[sa > 0].index
    tb = sb[sb > 0].index
    a2a = ta.intersection(tb).intersection(inds)
    a2a = ps[a2a].sort_values(inplace=False, ascending=True).index

    #BtoB
    ta = sa[sa < 0].index
    tb = sb[sb < 0].index
    b2b = ta.intersection(tb).intersection(inds)
    b2b = ps[b2b].sort_values(inplace=False, ascending=True).index
    return a2b, b2a, a2a, b2b


def plotChanges(data, na, nb, a2b, b2a, a2a, b2b, pcut, output):
    """
    Plot the compartment changes.
    """

    #plot the raw dots
    fig, ax = pylab.subplots(figsize=(3.2, 2.2))
    ax.scatter(data[na],
               data[nb],
               s=0.5,
               color="gray",
               alpha=0.6,
               label="total %s bins" % data.shape[0])

    #plot the changes
    ax.scatter(data[na][a2b],
               data[nb][a2b],
               s=1,
               color=colors[0],
               alpha=0.8,
               label="A->B %s bins" % len(a2b))
    ax.scatter(data[na][b2a],
               data[nb][b2a],
               s=1,
               color=colors[2],
               alpha=0.8,
               label="B->A %s bins" % len(b2a))
    ax.scatter(data[na][a2a],
               data[nb][a2a],
               s=1,
               color=colors[3],
               alpha=0.8,
               label="A->A %s bins" % len(a2a))
    ax.scatter(data[na][b2b],
               data[nb][b2b],
               s=1,
               color=colors[4],
               alpha=0.8,
               label="B->B %s bins" % len(b2b))

    leg = ax.legend(
        bbox_to_anchor=(1.05, 1.0),
        loc='upper left',
        labelcolor=["gray", colors[0], colors[2], colors[3], colors[4]])
    for h in leg.legendHandles:
        h._sizes = [10]
    ax.axvline(0, color="gray", linestyle="--")
    ax.axhline(0, color="gray", linestyle="--")
    ax.set_xlabel(f"{na} PC1")
    ax.set_ylabel(f"{nb} PC1")
    ax.set_title(f"Mahalanobis distance P-value < {pcut}")
    pylab.savefig(f"{output}_bins_flips.pdf")


def parseGtfFeature(t):
    ds = {}
    t = t.replace('"', '')
    for n in t.split('; '):
        s = n.split()
        ds[s[0]] = s[1]
    return ds


def readCompGenes(compf, gtf):
    """
    Read compartment and gene sets regions as HTSeq.GenomicArrayOfSets.
    """
    #read compartments
    comps = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for line in open(compf):
        line = line.split("\n")[0].split("\t")
        iv = HTSeq.GenomicInterval(line[0], int(line[1]), int(line[2]))
        comps[iv] += line[3]
    #read gtf
    genes = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gs = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        line = line.split("\n")[0].split("\t")
        if line[2] != 'exon':
            continue
        ds = parseGtfFeature(line[8])
        key = "|".join([ds["gene_id"], ds["gene_name"]])
        nline = [
            line[0], line[3], line[4],
            "|".join([ds["gene_id"], ds["gene_name"]]), ".", line[6]
        ]
        if key not in gs:
            gs[key] = [line[0], int(line[3]), int(line[4])]
        else:
            if int(line[3]) < gs[key][1]:
                gs[key][1] = int(line[3])
            if int(line[4]) > gs[key][2]:
                gs[key][2] = int(line[4])
    for g, v in gs.items():
        iv = HTSeq.GenomicInterval(v[0], v[1], v[2])
        genes[iv] += g
    return comps, genes


def anoBins(data, na, nb, s, comps, genes, fout):
    """
    Annotate changed bins.
    """
    ds = {}
    ags = set()
    for t in s:
        chrom = t.split(":")[0]
        start = t.split(":")[1].split("-")[0]
        end = t.split(":")[1].split("-")[1]
        iv = HTSeq.GenomicInterval(chrom, int(start), int(end))
        comp = list(list(comps[iv].steps())[0][1])[0]
        gs = set()
        for i, g in genes[iv].steps():
            gs.update(g)
        gs = list(gs)
        if len(gs) > 0:
            gs = [g.split("|")[-1] for g in gs]
        ags.update(gs)
        gs = ",".join(gs)
        ds[t] = {
            "chrom": chrom,
            "start": start,
            "end": end,
            f"{na} PC1": data.loc[t, na],
            f"{nb} PC1": data.loc[t, nb],
            "P-value": data.loc[t, "Chi-Square test P-value"],
            "compartmentId": comp,
            "overlappedGenes": gs,
        }
    ds = pd.DataFrame(ds).T
    ds.to_csv(fout + ".txt", index_label="binId", sep="\t")
    ags = list(ags)
    with open(fout + "_genes.list", "w") as fo:
        fo.write("\n".join(ags))


def compareComp(fa, fb, na, nb, gtf, output, pcut=0.01):
    """
    Pair-wise comparsion of compartments. 
    """
    print("Step 1: Stiching bins as compartments.")
    #step 1 get the compartment
    acompf = output + "_" + na + "_compartments.bed"
    bcompf = output + "_" + nb + "_compartments.bed"
    getCompartment(fa, acompf)
    getCompartment(fb, bcompf)

    print("Step 2: Performing two-passes mahalanobis distances caculation.")
    #step 2 prepare the bins level data
    sa = readBdg(fa)
    sb = readBdg(fb)
    t = sa.index.intersection(sb.index)
    sa = sa[t]
    sb = sb[t]
    data = pd.DataFrame({na: sa, nb: sb})

    #step 3 two passes MD test for bins
    dis, ps = twoPassesMDTest(data, pcut)
    inds = ps[ps < pcut].index

    #step 4 get flips or same compartment chanegs
    a2b, b2a, a2a, b2b = getFlips(sa, sb, ps, inds)

    print("Step 3: Plotting switch bins.")
    #step 5 show the changes
    plotChanges(data, na, nb, a2b, b2a, a2a, b2b, pcut, output)

    print("Step 4: Outputing results.")
    #step 6 output the p-values as bdg
    ps[ps < 1e-300] = 1e-300
    with open(f"{output}_-logP.bdg", "w") as fout:
        for i in ps.index:
            chrom = i.split(":")[0]
            start = i.split(":")[1].split("-")[0]
            end = i.split(":")[1].split("-")[1]
            line = [chrom, start, end, str(-np.log10(ps[i]))]
            fout.write("\t".join(line) + "\n")

    #step 7 output the distance and p-values
    data["Mahalanobis distance"] = dis
    data["Chi-Square test P-value"] = ps
    data.to_csv(f"{output}_MD_p-values.txt", sep="\t", index_label="binId")

    #step 8 annotate the changed bins and associated genes and output as bed files
    ds = {
        "AtoB": a2b,
        "BtoA": b2a,
        "AtoA": a2a,
        "BtoB": b2b,
    }
    comps, genes = readCompGenes(acompf, gtf)
    for k, v in ds.items():
        anoBins(data, na, nb, v, comps, genes, f"{output}_{k}_bins")


def main():
    op = help()
    compareComp(op.fa, op.fb, op.na, op.nb, op.gtf, op.output, op.pcut)


if __name__ == "__main__":
    main()
