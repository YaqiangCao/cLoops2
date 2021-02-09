import os,gzip
from glob import glob

import pandas as pd
from joblib import Parallel, delayed

#cLoops2 
from cLoops2.ds import PET
from cLoops2.utils import cFlush

def pre():
    fs = glob("../../13.distanceLinker/1.unqiue/*.bedpe.gz")
    ds = {}
    for f in fs:
        n = f.split("/")[-1].split(".bedpe")[0]
        fq = "../../11.filterFqs/%s_stat.txt"%n
        if not os.path.isfile(fq):
            print(f)
            continue
        peakBedpe = "../1.overlapWithPeaks/%s.bedpe.gz"%n
        ds[n] = {"bedpe":f,"fq":fq,"peakBedpe":peakBedpe}
    return ds


def get(bedpe,fq,peakbedpe,pre):
    #get reads in peaks
    readsInPeaks = set()
    for line in gzip.open(peakbedpe,"rt"):
        line = line.split("\n")[0].split("\t")
        readsInPeaks.add( line[6] )
    ds = pd.read_csv(fq,index_col=0,sep="\t")
    s = ds.sum(axis=1)
    t = s[s>-2]
    t = t.shape[0]/float(ds.shape[0])
    stat = {}
    for k in ["close","mid","distal","trans"]:
        stat[k] = {
                "inPeaksPETs": 0,
                "inPeaksLinker": 0,
                "notinPeaksPETs": 0, 
                "notinPeaksLinker": 0,
            }
    cis = 0
    for i,line in enumerate(gzip.open(bedpe,"rt")):
        if i % 100000 == 0:
            cFlush("%s PETs processed from %s" % (i, bedpe))
        line = line.split("\n")[0].split("\t")
        pet = PET(line)
        if s[line[6]] > -2:
            linker = True
        else:
            linker = False
        if line[6] in readsInPeaks:
            inPeak = True
        else:
            inPeak = False

        if pet.cis:
            cis += 1
            if pet.distance <= 1000:
                if inPeak:
                    stat["close"]["inPeaksPETs"] += 1
                    if linker:
                        stat["close"]["inPeaksLinker"] += 1
                else:
                    stat["close"]["notinPeaksPETs"] += 1
                    if linker:
                        stat["close"]["notinPeaksLinker"] += 1
            if 1000 < pet.distance <= 10000:
                if inPeak:
                    stat["mid"]["inPeaksPETs"] += 1
                    if linker:
                        stat["mid"]["inPeaksLinker"] += 1
                else:
                    stat["mid"]["notinPeaksPETs"] += 1
                    if linker:
                        stat["mid"]["notinPeaksLinker"] += 1
            if 10000 < pet.distance:
                if inPeak:
                    stat["distal"]["inPeaksPETs"] += 1
                    if linker:
                        stat["distal"]["inPeaksLinker"] += 1
                else:
                    stat["distal"]["notinPeaksPETs"] += 1
                    if linker:
                        stat["distal"]["notinPeaksLinker"] += 1
        else:
            if inPeak:
                stat["trans"]["inPeaksPETs"] += 1
                if linker:
                    stat["trans"]["inPeaksLinker"] += 1
            else:
                stat["trans"]["notinPeaksPETs"] += 1
                if linker:
                    stat["trans"]["notinPeaksLinker"] += 1
    rs = {  
            "totalRaw":ds.shape[0],
            "rawLinkerRatio(anyEnd)":t,
            "totalUnique":i,

            "transPETsRatio":(stat["trans"]["inPeaksPETs"]+stat["trans"]["notinPeaksPETs"])/float(ds.shape[0]),
            "transLinkerRatio": (stat["trans"]["inPeaksLinker"]+stat["trans"]["notinPeaksLinker"])/(stat["trans"]["inPeaksPETs"]+stat["trans"]["notinPeaksPETs"]),
            "transInPeaksLinkerRatio": stat["trans"]["inPeaksLinker"]/stat["trans"]["inPeaksPETs"],
            "transNotInPeaksLinkerRatio": stat["trans"]["notinPeaksLinker"]/stat["trans"]["notinPeaksPETs"],

            "closePETsRatio(inCis)":(stat["close"]["inPeaksPETs"]+stat["close"]["notinPeaksPETs"])/cis,
            "closeLinkerRatio": (stat["close"]["inPeaksLinker"]+stat["close"]["notinPeaksLinker"])/(stat["close"]["inPeaksPETs"]+stat["close"]["notinPeaksPETs"]),
            "closeInPeaksLinkerRatio": stat["close"]["inPeaksLinker"]/stat["close"]["inPeaksPETs"],
            "closeNotInPeaksLinkerRatio": stat["close"]["notinPeaksLinker"]/stat["close"]["notinPeaksPETs"],

            "midPETsRatio":(stat["mid"]["inPeaksPETs"]+stat["mid"]["notinPeaksPETs"])/cis,
            "midLinkerRatio": (stat["mid"]["inPeaksLinker"]+stat["mid"]["notinPeaksLinker"])/(stat["mid"]["inPeaksPETs"]+stat["mid"]["notinPeaksPETs"]),
            "midInPeaksLinkerRatio": stat["mid"]["inPeaksLinker"]/stat["mid"]["inPeaksPETs"],
            "midNotInPeaksLinkerRatio": stat["mid"]["notinPeaksLinker"]/stat["mid"]["notinPeaksPETs"],
            
            "distalPETsRatio":(stat["distal"]["inPeaksPETs"]+stat["distal"]["notinPeaksPETs"])/cis,
            "distalLinkerRatio": (stat["distal"]["inPeaksLinker"]+stat["distal"]["notinPeaksLinker"])/(stat["distal"]["inPeaksPETs"]+stat["distal"]["notinPeaksPETs"]),
            "distalInPeaksLinkerRatio": stat["distal"]["inPeaksLinker"]/stat["distal"]["inPeaksPETs"],
            "distalNotInPeaksLinkerRatio": stat["distal"]["notinPeaksLinker"]/stat["distal"]["notinPeaksPETs"],

    }
    rs = pd.Series(rs)
    print(pre)
    print(rs)
    return pre,rs


ds = pre()
ds = Parallel(n_jobs=20)(delayed(get)(v["bedpe"],v["fq"],v["peakBedpe"],k) for k, v in ds.items())
data = {}
for d in ds:
    data[d[0]] = d[1]
data = pd.DataFrame(data).T
data = data[ [ "totalRaw", "rawLinkerRatio(anyEnd)", "totalUnique", "transPETsRatio", "transLinkerRatio", "transInPeaksLinkerRatio", "transNotInPeaksLinkerRatio", "closePETsRatio(inCis)", "closeLinkerRatio", "closeInPeaksLinkerRatio", "closeNotInPeaksLinkerRatio", "midPETsRatio", "midLinkerRatio", "midInPeaksLinkerRatio", "midNotInPeaksLinkerRatio", "distalPETsRatio", "distalLinkerRatio", "distalInPeaksLinkerRatio", "distalNotInPeaksLinkerRatio", ] ]
data.to_csv("linker_stat_inpeaks.txt",sep="\t")
