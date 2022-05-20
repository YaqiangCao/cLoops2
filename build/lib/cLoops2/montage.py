#!/usr/bin/env python
#--coding:utf-8 --
"""
cLoops2 montage module. 
Mostly inspired by https://github.com/tsznxx/PyCircos. 
2020-08-10: bascically finished calss:Circos, may need updates for each PETs as arc
2020-08-11: updated as line width for strength
2020-08-12: try to highlight target specific regions, basically finished. If following time available, add text annotation. 1D data bin/smooth needed. 
2020-08-16: updated as adding single region plot
2021-04-02: add option no 1D signal
2021-04-28: update view point mode
"""
#sys
import json

#3rd
import numpy as np
import pandas as pd
from matplotlib.path import Path
from matplotlib.patches import PathPatch

#cLoops2
from cLoops2.ds import XY
from cLoops2.io import parseIxy, parseTxt2Loops
from cLoops2.cmat import get1DSig
from cLoops2.settings import *


class Circos(object):
    '''
    Python Circos.
    '''

    def __init__(self, regions, length='length', figsize=(4, 4), gap=0.5):
        '''
        Initiation from a list of beds.
        Parameters:
            regions: region sizes as pd.Dataframe
            size: string
                column name for chrom sizes
            figsize: tuple or list of floats
                width x height
            gap: float
                gap between each Bed.  
        '''
        self.regions = regions
        total_len = self.regions["length"].sum()
        self.len_per_degree = total_len / (360. - gap * self.regions.shape[0])
        self.len_per_theta = total_len / (
            np.pi * 2 - np.deg2rad(gap) * self.regions.shape[0])
        # chromosome start thetas
        cumlen = [0] + list(
            self.regions.length.cumsum())[:-1]  # accumulative length
        self.regions['theta_start'] = [
            np.deg2rad(l / self.len_per_degree + gap * i)
            for i, l in enumerate(cumlen)
        ]
        # polar axis
        self.fig = pylab.figure(figsize=figsize)
        l = max(figsize)
        self.pax = self.fig.add_axes([0, 0, 1, 1], polar=True)
        self.pax.axis('off')
        self.pax.set_ylim(0, l)

    def get_theta(self, gid, pos):
        '''
        get the theta of the position.
        Parameters:
            gid: string, chrom labels
            pos: int, chrom coordinates
        '''
        et = self.regions.loc[gid, 'theta_start'] + pos / self.len_per_theta
        return et

    def draw_scaffold(self, rad, width, colors=[], fill=False, **kwargs):
        '''
        Draw scaffold.
        Parameters:
            rad: float
                radius.
            width: float
                width of the band. negative width means inner width. eg. rad=8,width=1 equal to rad=9,width=-1.
            colors: list of colors
                cycling colors. at least two colors.
            alpha: float
                alpha value.
        '''
        n = len(colors)
        if fill == False or n == 0:
            kwargs.update({
                'edgecolor': 'k',
                'linewidth': 1,
                'linestyle': '-',
                'fill': False
            })
        else:
            kwargs.update({'linewidth': 0})
        for i, gid in enumerate(self.regions.index):
            if n:
                kwargs['color'] = colors[i]
            et1, et2 = self.regions.theta_start[gid], self.get_theta(
                gid, self.regions.length[gid] - 1)
            self.pax.bar([(et1 + et2) / 2], [width],
                         width=et2 - et1,
                         bottom=rad,
                         **kwargs)

    def draw_scaffold_ids(self, rad, text=None, inside=False, **kwargs):
        '''
        Draw scaffold region IDs.
        Parameters:
            rad: float
                radius
            inside: bool
                draw chrom labels inside
            kwargs: dict
                to ax.annotate()
                    fontsize, rotation
        '''
        kwargs.setdefault('ha', 'center')
        kwargs.setdefault('va', 'center')
        rotation = kwargs.get('rotation', 0)
        ml = max([len(gid) for gid in self.regions.index])
        for gid in self.regions.index:
            deg = np.rad2deg(self.get_theta(gid, self.regions.length[gid] / 2))
            kwargs['rotation'] = rotation + deg
            if 90 < kwargs['rotation'] < 270:
                kwargs['rotation'] += 180
            if inside:  # add spaces to the right side
                lstr = ' ' * (ml - len(gid)) + gid
            else:
                lstr = gid + ' ' * (ml - len(gid))
            if text is not None:
                t = gid + " " + str(text[gid])
            else:
                t = gid
            self.pax.annotate(t, xy=[np.deg2rad(deg), rad], **kwargs)

    def draw_link(self,
                  rad,
                  gids,
                  start,
                  end,
                  color=None,
                  label=None,
                  alpha=1,
                  lw=1):
        '''
        Draw links
        Parameters:
            rad: float,radius
            gids: list,list of two chroms
            starts, ends: list,list of start/end coordinates
            color: string,face color
            alpha: float alpha            
        '''
        ets = self.get_theta(gids[0], start)
        ete = self.get_theta(gids[1], end)
        points = [
            (ets, rad),
            (0, 0),
            (ete, rad),
        ]
        # parse patches
        codes = [Path.CURVE3] * len(points)
        codes[0] = Path.MOVETO
        path = Path(points, codes)
        if color is None:
            color = "k"
        patch = PathPatch(path,
                          alpha=alpha,
                          lw=lw,
                          facecolor="none",
                          edgecolor=color,
                          label=label)
        self.pax.add_patch(patch)

    def fill_between(self,
                     rad,
                     data,
                     scale=1,
                     gid='chrom',
                     start='start',
                     end='end',
                     score='score',
                     color="red",
                     vmin=None,
                     vmax=None,
                     alpha=1,
                     label=""):
        '''
        Draw densities.
        Parameters:
            rad: float, radius
            data: pandas.DataFrame,chromosomal regions
            start, end: int, chrom start or end
            score: float,chrom interval scores
            cutoff: float,abs(value) < cutoff are filled in grey
            scale: float,scalling factor of original scores
        The data should be normalized first
        '''
        if vmax is None:
            vmax = data[score].max()
        else:
            vmax = vmax
        if vmin is None:
            vmin = 0
        else:
            vmin = vmin
        height = scale / (vmax - vmin)
        i = 0
        for gid, start, end, score in zip(data[gid], data[start], data[end],
                                          data[score]):
            ets = self.get_theta(gid, start)
            ete = self.get_theta(gid, end)
            if score <= vmin:
                continue
            if score >= vmax:
                score = vmax
            h = (score - vmin) * height
            if i == 0:
                self.pax.fill_between([ets, ete], [rad, rad],
                                      [rad + h, rad + h],
                                      color=color,
                                      label=label,
                                      alpha=alpha)
                i = i + 1
            else:
                self.pax.fill_between([ets, ete], [rad, rad], [rad + h],
                                      rad + h,
                                      color=color,
                                      alpha=alpha)


def mergeCov(cov):
    ncov = {}
    k = 0
    i = 0
    while i < len(cov) - 1:
        if cov[i] == 0:  #find the non 0 start
            i += 1
            continue
        for j in range(i + 1, len(cov)):  #find the same value stop
            if cov[j] != cov[i]:
                break
        ncov[k] = {"start": i, "end": j, "cov": cov[i]}
        k = k + 1
        if j == len(cov) - 1:
            break
        i = j
    return ncov


def montage(
        fixy,
        bed,
        fout,
        ext=5,
        viewPoint="",
        cut=0,
        mcut=-1,
        vmin=None,
        vmax=None,
        simple=True,
        ppmw=10,
        aw=0.5,
        oneD=True,
):
    """
    Montage analysis of specific regions
    @param fixy: string, .ixy file path
    @param bed: string, .bed file path, 4th column should be region names
    @param fout: string, output file prefix
    @param ext: int, extesion fold for up-stream and down-stream 
    @param vmin: float, minal scale for 1D 
    @param vmax: float, maxial scale for 1D
    @param ppmw: int, 1 PETs per million width
    @param oneD: bool, whether to plot 1D profile. For data like Hi-C, better not.
    """
    #step 0, basic settings
    rad = 2.5
    gap = 3
    scaffoldw = 0.15
    oneDscale = 1
    if viewPoint != "":
        if "," in viewPoint:
            viewPoint = viewPoint.split(",")
        else:
            viewPoint = [viewPoint]
    #step 1, parse regions
    regions = {}
    for line in open(bed):
        line = line.split("\n")[0].split("\t")
        rid = line[3]
        s = int(line[1])
        e = int(line[2])
        w = e - s
        regions[rid] = {
            "rawStart": s,
            "rawEnd": e,
            "extStart": s - ext * w,
            "extEnd": e + ext * w,
            "newStart": 0,
            "newEnd": (2 * ext + 1) * w,
            "centerStart": ext * w,
            "centerEnd": (ext + 1) * w,
            "length": (2 * ext + 1) * w,
        }
    if len(regions) == 0 :
        print("No regions in input bed file, return")
        return
    regions = pd.DataFrame(regions).T
    if regions.shape[0] == 1:
        print("Only 1 region in input bed file, all interactions in the region will be ploted.") 

    #step 2, plot basic scaffold
    cr = Circos(regions, gap=gap)
    cr.draw_scaffold(rad,
                     scaffoldw,
                     colors=["gray"] * regions.shape[0],
                     fill=True)
    #highlight the target region
    for i, rid in enumerate(regions.index):
        thetas = cr.get_theta(rid, regions.loc[rid, "centerStart"])
        thetae = cr.get_theta(rid, regions.loc[rid, "centerEnd"])
        theta = (thetas + thetae) / 2
        width = thetae - thetas
        if i == 0:
            cr.pax.bar(theta,
                       scaffoldw,
                       width=width,
                       bottom=rad,
                       color=colors[0],
                       alpha=0.8,
                       label="target region, up/down-stream extend %s fold" %
                       ext)
        else:
            cr.pax.bar(theta,
                       scaffoldw,
                       width=width,
                       bottom=rad,
                       color=colors[0],
                       alpha=0.8)

    #step 3, plot interactions as arc
    metaf = "/".join(fixy.split("/")[:-1]) + "/petMeta.json"
    meta = json.loads(open(metaf).read())
    #total reads as million
    tot = meta["Unique PETs"] / 10**6
    chrom, xy = parseIxy(fixy, cut=cut, mcut=mcut)
    xy2 = XY(xy[:, 0], xy[:, 1])  #XY object
    
    #step 3, if only 1 regions, plot all archs
    if regions.shape[0] == 1:
        ra = regions.index[0]
        cs = xy2.queryPeakBoth(regions.loc[ra]["extStart"],
                                regions.loc[ra]["extEnd"])
        cs = list(cs)
        for x, y in xy[cs, ]:
            nx = x - regions.loc[ra, "extStart"]
            ny = y - regions.loc[ra, "extStart"]
            cr.draw_link(rad, [ra, ra],
                         nx,
                         ny,
                         lw=0.5,
                         alpha=0.5,
                         color="purple")
        cr.draw_link(0, 
                    [ra, ra],
                    0,
                    0,
                    lw=0.5,
                    alpha=0.5,
                    color="purple",
                    label="1 PET")
 
    else:
    #step 4.1, plot representative
        if simple:
            data = {}
            for ra in regions.index:
                data[ra] = {}
                for rb in regions.index:
                    if rb not in data:
                        data[rb] = {}
                    if ra == rb:
                        continue
                    ca, cb, cab = xy2.queryLoop(regions.loc[ra]["rawStart"],
                                                regions.loc[ra]["rawEnd"],
                                                regions.loc[rb]["rawStart"],
                                                regions.loc[rb]["rawEnd"])
                    data[ra][rb] = len(cab)
                    data[rb][ra] = len(cab)
            data = pd.DataFrame(data)
            data = data.fillna(0)
            ns = list(data.index)
            ns.sort()
            data = data.loc[ns, ns]
            data = data.astype("float")
            data = data / tot
            data.to_csv("%s_interactionPETsPerMillion.txt" % fout, sep="\t")
            for na in data.index:
                if viewPoint != "" and na not in viewPoint:
                    continue
                for nb in data.columns:
                    if data.loc[na, nb] == 0.0:
                        continue
                    lw = data.loc[na, nb] * ppmw
                    s = (regions.loc[na, "centerStart"] +
                         regions.loc[na, "centerEnd"]) / 2
                    e = (regions.loc[nb, "centerStart"] +
                         regions.loc[nb, "centerEnd"]) / 2
                    cr.draw_link(rad, [na, nb],
                                 s,
                                 e,
                                 lw=lw,
                                 alpha=0.8,
                                 color="purple")
            cr.draw_link(0, [na, na],
                         0,
                         0,
                         #lw=1,
                         #label="%s PETs per million"%ppmw,
                         lw=ppmw,
                         label="1 PETs per million",
                         alpha=0.8,
                         color="purple")
        #step 4.2, plot all archs
        else:
            cors = set()
            for ra in regions.index:
                if viewPoint != "" and ra not in viewPoint:
                    continue
                for rb in regions.index:
                    if ra == rb:
                        continue
                    if viewPoint == "":
                        ### important!!! always keep lefta,righta, leftb,rightb, righta <=leftb
                        if regions.loc[rb, "rawStart"] < regions.loc[ra, "rawEnd"]:
                            continue
                        ca, cb, cab = xy2.queryLoop(regions.loc[ra,"rawStart"],
                                                    regions.loc[ra,"rawEnd"],
                                                    regions.loc[rb,"rawStart"],
                                                    regions.loc[rb,"rawEnd"])
                        cab = list(cab)
                        for x, y in xy[cab, ]:
                            nx = x - regions.loc[ra, "extStart"]
                            ny = y - regions.loc[rb, "extStart"]
                            cr.draw_link(rad, [ra, rb],
                                         nx,
                                         ny,
                                         #lw=0.5,
                                         lw=aw,
                                         alpha=0.5,
                                         color="purple")
                    else:
                        ca, cb, cab = xy2.queryLoop(regions.loc[ra,"rawStart"],
                                                    regions.loc[ra,"rawEnd"],
                                                    regions.loc[rb,"rawStart"],
                                                    regions.loc[rb,"rawEnd"])
                        cab = list(cab)
                        cab = [ t for t in cab if t not in cors]
                        cors.update(cab)
                        if regions.loc[ra,"rawStart"] < regions.loc[rb,"rawStart"]:
                            for x, y in xy[cab, ]:
                                nx = x - regions.loc[ra, "extStart"]
                                ny = y - regions.loc[rb, "extStart"]
                                cr.draw_link(rad, [ra, rb],
                                             nx,
                                             ny,
                                             #lw=0.5,
                                             lw=aw,
                                             alpha=0.5,
                                             color="purple")
                        else:
                            for x, y in xy[cab, ]:
                                nx = x - regions.loc[rb, "extStart"]
                                ny = y - regions.loc[ra, "extStart"]
                                cr.draw_link(rad, [rb, ra],
                                             nx,
                                             ny,
                                             #lw=0.5,
                                             lw=aw,
                                             alpha=0.5,
                                             color="purple")
            cr.draw_link(0, [ra, ra],
                         0,
                         0,
                         lw=0.5,
                         alpha=0.5,
                         color="purple",
                         label="1 PET")
    
    if oneD:
        #step 5, plot 1D
        covData = {}
        i = 0
        for rid in regions.index:
            s = get1DSig(xy2, int(regions.loc[rid, "extStart"]),
                         int(regions.loc[rid, "extEnd"]))
            s = mergeCov(s)
            for k, v in s.items():
                covData[i] = {
                    "rid": rid,
                    "start": v["start"],
                    "end": v["end"],
                    "cov": v["cov"]
                }
                i += 1
        #normalization the 1D signal
        covData = pd.DataFrame(covData).T
        s = covData["cov"] / tot
        covData["cov"] = s
        if vmin is not None:
            vmin = vmin
        else:
            vmin = 0
        if vmax is not None:
            vmax = vmax
        else:
            vmax = covData["cov"].max()
        label = "1D scale (RPM):[%.3f, %.3f]" % (vmin, vmax)
        cr.fill_between(rad + scaffoldw,
                        covData,
                        gid="rid",
                        start="start",
                        end="end",
                        score="cov",
                        color=colors[1],
                        alpha=0.8,
                        scale=oneDscale,
                        vmin=vmin,
                        vmax=vmax,
                        label=label)
        rad = rad + scaffoldw + oneDscale
    else:
        rad = rad + scaffoldw

    #step 6, plot region ids
    cr.draw_scaffold_ids(rad, inside=False, fontsize=8)
    #save the legends and adjust
    leg = cr.pax.legend(bbox_to_anchor=(0.8, -0.05))
    #line0 = list(leg.get_lines())[0]
    #pylab.setp(line0,width=ppmw)
    #save the results
    pylab.savefig('%s_rehoboam.pdf' % fout)
