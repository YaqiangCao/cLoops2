#!/usr/bin/env python3
#--coding:utf-8 --
"""
ano.py 

Previouse analyzeLoops.py, now changed to ano.py to also annotate peaks/domains. 

Include cLoops2 loops-centric analysis module. Input should be the _loops.txt file and other annotation files. Mainly contain following analysis. 
- [x] loops annotation to target genes as enhancer and promoter 
- [x] loops annotation to target genes through network method
- [x] find HUBs through HITS algorithm

2020-11-04: modified to first find overlapped TSS, if no or multiple, then find the closest one.
2021-01-21: peaks annotation going to be added.
"""

#sys
import os

#3rd
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from scipy.spatial import KDTree
from joblib import Parallel, delayed

#cLoops2
from cLoops2.ds import Exon, Gene, Peak
from cLoops2.io import parseTxt2Loops


def parseGtfLine(line, tid=False):
    """
    Parse gene gtf line.
    """
    e = Exon()
    e.chrom = line[0]
    e.start = int(line[3])
    e.end = int(line[4])
    e.length = e.end - e.start
    e.strand = line[6]
    attr = line[8].replace('"', '').split(";")
    ts = {}
    for t in attr:
        t = t.split()
        if len(t) != 2:
            continue
        ts[t[0]] = t[1]
    if tid:
        try:
            e.name = ts["transcript_name"]
        except:
            e.name = ts["gene_name"]
    else:
        e.name = ts["gene_name"]
    if tid:
        e.id = ts["transcript_id"]
    else:
        e.id = ts["gene_id"]
    return e


def readGenes(gtf, tid=False):
    """
    Read gene annotion file for genes or transcripts
    """
    gs = {}
    #get all genes information
    print("reading annotaions from %s" % gtf)
    for line in tqdm(open(gtf).read().split("\n")[:-1]):
        if line.startswith("#"):
            continue
        line = line.split("\n")[0].split("\t")
        if line[2] != "exon":
            continue
        e = parseGtfLine(line, tid)
        if e.name not in gs:
            g = Gene()
            g.chrom = e.chrom
            g.start = e.start
            g.end = e.end
            g.strand = e.strand
            g.name = e.name
            g.id = e.id
            g.exons = {(e.start, e.end): e}
            gs[g.name] = g
        else:
            #same position exons
            if (e.start, e.end) in gs[e.name].exons:
                continue
            else:
                g = gs[e.name]
                if e.start < g.start:
                    g.start = e.start
                if e.end > g.end:
                    g.end = e.end
                g.exons[(e.start, e.end)] = e
    #get all genes information
    ngs = {}  #key is chromosome
    for k, g in gs.items():
        if g.chrom not in ngs:
            ngs[g.chrom] = {}
        if g.strand == "+":
            tss = g.start
        else:
            tss = g.end
        #tss position is key, other information is value, for following search
        if tss not in ngs[g.chrom]:
            ngs[g.chrom][tss] = g
    return ngs


def findOverlapOrNearest(gs, ts, tree, start, end):
    """
    first to check direct overlap with TSS, if no or multiple items, then get the close one
    @param gs: {tss: cLoops2.ds.Gene}, tss is key and int
    @pram ts: [tss]
    @param tree: KDTree from TSSs
    @param start: query start
    @param end: query end
    return gene and distance
    """
    #step 1, find overlaps
    rs = set()
    for i in range(start, end + 1):
        if i in gs:
            rs.add(gs[i])
    if len(rs) > 0:
        rs = list(rs)
        return rs, [0] * len(rs)
    #find the nearest one
    else:
        d, i = tree.query([(start + end) / 2], k=1)
        g = gs[ts[i][0]]
        #d = ts[i][0] - (start+end)/2
        d = int(d)
        return [g], [d]


def findNearestTss(chrom, loops, gs, pdis=2000):
    """
    Find nearest TSS for loop anchors. 
    @param loops: list of cLoops2.ds.Loop
    @param gs: {tss: cLoops2.ds.Gene}, tss is key and int
    """
    ts = np.array([[tss] for tss in gs.keys()])
    cov = {}
    for tss, g in gs.items():
        cov[tss] = g
    tree = KDTree(ts)
    ds = {}
    for loop in loops:
        xgs, xds = findOverlapOrNearest(gs, ts, tree, loop.x_start, loop.x_end)
        ygs, yds = findOverlapOrNearest(gs, ts, tree, loop.y_start, loop.y_end)
        if len(xgs) > 1:
            xt = "Promoter"
            xd = 0
        else:
            xd = xds[0]
            if abs(xd) <= pdis:
                xt = "Promoter"
            else:
                xt = "Enhancer"
        if len(ygs) > 1:
            yt = "Promoter"
            yd = 0
        else:
            yd = yds[0]
            if abs(yd) <= pdis:
                yt = "Promoter"
            else:
                yt = "Enhancer"
        ds[loop.id] = {
            "typeAnchorA":
            xt,
            "typeAnchorB":
            yt,
            "nearestDistanceToGeneAnchorA":
            xd,
            "nearestDistanceToGeneAnchorB":
            yd,
            "nearestTargetGeneAnchorA":
            ",".join([
                xg.chrom + ":" + str(xg.start) + "-" + str(xg.end) + "|" +
                xg.strand + "|" + xg.name for xg in xgs
            ]),
            "nearestTargetGeneAnchorB":
            ",".join([
                yg.chrom + ":" + str(yg.start) + "-" + str(yg.end) + "|" +
                yg.strand + "|" + yg.name for yg in ygs
            ]),
        }
    return ds


def annotateLoopToGenes(loops, genes, fout, pdis=2000, cpu=1):
    """
    Annotate loops releative to genes. 
    @param loops: { "chrom-chrom":[] }, in list are cLoops2.ds.Loop
    """
    print("Annotating loops to enhancers and promoters.")
    ks = [key for key in loops.keys() if key in genes]
    ds = Parallel(n_jobs=cpu,
                  backend="multiprocessing")(delayed(findNearestTss)(
                      chrom,
                      loops[chrom],
                      genes[chrom],
                      pdis=pdis,
                  ) for chrom in tqdm(ks))
    rs = {}
    for d in ds:
        for k, v in d.items():
            rs[k] = v
    rs = pd.DataFrame(rs).T
    fo = fout + "_LoopsGtfAno.txt"
    rs.to_csv(fo, sep="\t", index_label="loopId")
    return fo


def stichAnchors(chrom, loops, margin=1):
    """
    Stich close anchors based on postion array.
    """
    cov = set()
    for i, loop in enumerate(loops):
        cov.update(range(loop.x_start, loop.x_end + 1))
        cov.update(range(loop.y_start, loop.y_end + 1))
    cov = list(cov)
    cov.sort()
    npeaks = []
    i = 0
    while i < len(cov) - 1:
        j = i + 1
        while j < len(cov):
            if cov[j] - cov[j - 1] > margin:
                break
            else:
                j += 1
        peak = Peak()
        peak.chrom = chrom
        peak.start = cov[i]
        peak.end = cov[j - 1]
        peak.length = cov[j - 1] - cov[i] + 1
        npeaks.append(peak)
        i = j  #update search start
    return npeaks


def getNet(chrom, loops, genes, pdis=2000, gap=1):
    """
    Get the enhancer/promoter network for one chromosome.  
    """
    #step 1 get merged anchors
    anchors = stichAnchors(chrom, loops, margin=gap)
    #step 2 annotate anchors
    nanchors = {}
    ts = np.array([[tss] for tss in genes.keys()])
    tree = KDTree(ts)
    for anchor in anchors:
        gs, ds = findOverlapOrNearest(genes, ts, tree, anchor.start,
                                      anchor.end)
        if len(gs) > 1:
            t = "Promoter"
            d = 0
        else:
            d = ds[0]
            if abs(d) <= pdis:
                t = "Promoter"
            else:
                t = "Enhancer"
        n = anchor.chrom + ":" + str(anchor.start) + "-" + str(
            anchor.end) + "|" + t
        nanchors[n] = {
            "chrom":
            anchor.chrom,
            "start":
            anchor.start,
            "end":
            anchor.end,
            "type":
            n.split("|")[-1],
            "nearestDistanceToTSS":
            d,
            "nearestGene":
            ",".join([g.name for g in gs]),
            "nearestGeneLoc":
            ",".join([
                g.chrom + ":" + str(g.start) + "-" + str(g.end) + "|" +
                g.strand + "|" + g.name for g in gs
            ])
        }
    anchors = nanchors
    del nanchors
    #step 3 assign each anchor to merged annotated anchor and build the network
    anchorCov = {}
    for k, v in anchors.items():
        for i in range(v["start"], v["end"] + 1):
            anchorCov[i] = k
    ds = {}  #anchor annotations
    nets = {}  #net information
    G = nx.Graph()  #networkx graph structure
    for loop in loops:
        xt, yt = None, None
        for i in range(loop.x_start, loop.x_end + 1):
            if i in anchorCov:
                xt = anchorCov[i]
                break
        for i in range(loop.y_start, loop.y_end + 1):
            if i in anchorCov:
                yt = anchorCov[i]
                break
        ds[loop.id] = {
            "mergedAnchorA": xt,
            "mergedAnchorB": yt,
        }
        if xt == yt:
            continue
        ns = [xt, yt]
        ns.sort()  #sort for converging keys
        if ns[0] not in nets:
            nets[ns[0]] = set()
        nets[ns[0]].add(ns[1])
        #network edges
        G.add_edge(ns[0], ns[1])
    #step 4 find all enhancers linked to target gene
    targets = {}
    #step 4.1 find the direct enhancer that link to promoter
    for node in G.nodes:
        if node.split("|")[-1] == "Promoter":
            if node in targets:
                continue
            targets[node] = {
                "targetGene": anchors[node]["nearestGeneLoc"],
                "directEnhancer": set(),
                "indirectEnhancer": set(),
                "directPromoter": set(),
                "indirectPromoter": set()
            }
            ns = list(nx.descendants(G, node))
            #find all releated nodes
            for n in ns:
                p = nx.algorithms.shortest_path(G, source=node, target=n)
                if n.split("|")[-1] == "Promoter":
                    if len(p) == 2:
                        targets[node]["directPromoter"].add(n)
                    else:
                        targets[node]["indirectPromoter"].add(n)
                if n.split("|")[-1] == "Enhancer":
                    if len(p) == 2:
                        targets[node]["directEnhancer"].add(n)
                    else:
                        targets[node]["indirectEnhancer"].add(n)
            #step 4.2. find hub enhancer
            #only using non-redundant node to find hubs
            nns = []
            tmp = set()
            for n in ns:
                tn = n.split("|")[0]
                if tn not in tmp:
                    nns.append(n)
                tmp.add(tn)
            ns = list(nns)
            ns.append(node)
            subg = G.subgraph(ns)
            try:
                hubs, authorities = nx.hits(subg,
                                            max_iter=1000,
                                            normalized=True)
            except:
                print(
                    "For %s, hard to find the hub by running HITS algorithm of 1000 iteration."
                    % node)
                targets[node]["directEnhancerHub"] = ""
                targets[node]["indirectEnhancerHub"] = ""
                continue
            hubs = pd.Series(hubs)
            hubs = hubs.sort_values(inplace=False, ascending=False)
            if len(targets[node]["directEnhancer"]) >= 2:
                des = hubs[list(targets[node]["directEnhancer"])]
                des = des.sort_values(inplace=False, ascending=False)
                targets[node]["directEnhancerHub"] = des.index[0]
            else:
                targets[node]["directEnhancerHub"] = ""
            if len(targets[node]["indirectEnhancer"]) >= 2:
                indes = hubs[list(targets[node]["indirectEnhancer"])]
                indes = indes.sort_values(inplace=False, ascending=False)
                targets[node]["indirectEnhancerHub"] = indes.index[0]
            else:
                targets[node]["indirectEnhancerHub"] = ""
    return anchors, ds, nets, targets


def getNetworksFromLoops(loops, genes, fout, pdis=2000, gap=1, cpu=1):
    """
    Merge overlapped acnhors first then construct interaction network.
    """
    ks = [key for key in loops.keys() if key in genes]
    print("Merging anchors and annotating loops through networks.")
    ds = Parallel(n_jobs=cpu, backend="multiprocessing")(delayed(getNet)(
        chrom,
        loops[chrom],
        genes[chrom],
        pdis=pdis,
        gap=gap,
    ) for chrom in tqdm(ks))
    anchors, anots, nets, targets = {}, {}, {}, {}
    for d in ds:
        for k, v in d[0].items():
            anchors[k] = v
        for k, v in d[1].items():
            anots[k] = v
        for k, v in d[2].items():
            nets[k] = v
        for k, v in d[3].items():
            targets[k] = v
    #output results
    #anchors
    anchors = pd.DataFrame(anchors).T
    anchors.to_csv(fout + "_mergedAnchors.txt", sep="\t", index_label="anchor")
    with open(fout + "_mergedAnchors.bed", "w") as fo:
        for t in anchors.itertuples():
            line = [t[1], t[2], t[3], t[0]]
            fo.write("\t".join(list(map(str, line))) + "\n")
    #annotations
    anots = pd.DataFrame(anots).T
    anots.to_csv(fout + "_loop2anchors.txt", sep="\t", index_label="loopId")
    #networks
    with open(fout + "_ep_net.sif", "w") as fo:
        for s, es in nets.items():
            es = list(es)
            ta = s.split("|")[-1]
            for e in es:
                tb = e.split("|")[-1]
                t = [ta, tb]
                t.sort()
                t = "-".join(t)
                line = [s, t, e]
                fo.write("\t".join(line) + "\n")
    with open(fout + "_targets.txt", "w") as fo:
        ks = list(targets.keys())
        ks.sort()
        line = [
            "Promoter", "PromoterTarget", "directEnhancer", "indirectEnhancer",
            "directPromoter", "indirectPromoter", "directEnhancerHub",
            "indirectEnhancerHub"
        ]
        fo.write("\t".join(line) + "\n")
        for k in ks:
            line = [
                k, targets[k]["targetGene"],
                ",".join(targets[k]["directEnhancer"]),
                ",".join(targets[k]["indirectEnhancer"]),
                ",".join(targets[k]["directPromoter"]),
                ",".join(targets[k]["indirectPromoter"]),
                targets[k]["directEnhancerHub"],
                targets[k]["indirectEnhancerHub"]
            ]
            fo.write("\t".join(line) + "\n")

### annotate loops
def anaLoops(loopf,
             fout,
             gtf=None,
             tid=False,
             pdis=2000,
             net=False,
             gap=1,
             cpu=1):
    """
    Analyze loops.
    @param loopf: str, name of loops file,  _loops.txt or _dloops.txt file
    @param fout: str, output prefix
    @param gtf: str, GTF file name 
    @param tid: bool, if set true, use transcript id for alternative TSS
    @param pdis: <=distance nearest TSS to define as promoter
    @param net: bool, whether use network search for all linked anchors/enhancers/promoters for target gene
    @param gap: int, gap for merge anchors
    @param cpu: int, number of CPU to run analysis
    """
    loops = parseTxt2Loops(loopf, cut=0)
    #only annotate cis loops
    nloops = {}
    for key in loops.keys():
        nk = key.split("-")
        if nk[0] != nk[1]:
            continue
        nloops[nk[0]] = loops[key]
    loops = nloops
    if gtf is not None and gtf != "":
        if not os.path.isfile(gtf):
            print("Input %s not exists, continue to other analysis." % gtf)
        else:
            #gene annotions, {chrom:{tss:g}}, tss is int
            genes = readGenes(gtf, tid=tid)
            anf = annotateLoopToGenes(loops, genes, fout, pdis=pdis, cpu=cpu)
            #get common summary of interaction type summary and distance summary
            if net:
                getNetworksFromLoops(loops,
                                     genes,
                                     fout,
                                     pdis=pdis,
                                     gap=gap,
                                     cpu=cpu)
### annotate peaks
def anaPeaks(peakf,
             fout,
             gtf=None,
             tid=False,
             pdis=2000,
             gap=1,
             cpu=1):
    """
    Annotate peaks.
    @param peakf: str, name of loops file,  _loops.txt or _dloops.txt file
    @param fout: str, output prefix
    @param gtf: str, GTF file name 
    @param tid: bool, if set true, use transcript id for alternative TSS
    @param pdis: <=distance nearest TSS to define as promoter
    @param net: bool, whether use network search for all linked anchors/enhancers/promoters for target gene
    @param gap: int, gap for merge anchors
    @param cpu: int, number of CPU to run analysis
    """
    if gtf is not None and gtf != "":
        if not os.path.isfile(gtf):
            print("Input %s not exists, continue to other analysis." % gtf)
        else:
            #gene annotions, {chrom:{tss:g}}, tss is int
            genes = readGenes(gtf, tid=tid)
