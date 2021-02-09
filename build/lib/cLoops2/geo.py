#!/usr/bin/env python
#--coding:utf-8--
"""
geo.py
Geometric methods for intervals (anchors) and rectangles (loops).
"""

from copy import deepcopy
from cLoops2.ds import Peak


#peaks related operations
def checkPeakOverlap(peaka, peakb, margin=None):
    """
    check the overlap of a region for the same chromosome
    """
    if peaka.chrom != peakb.chrom:
        return False
    if peakb.start <= peaka.start <= peakb.end or peakb.start <= peaka.end <= peakb.end:
        return True
    if peaka.start <= peakb.start <= peaka.end or peaka.start <= peakb.end <= peaka.end:
        return True
    if margin is not None:
        if peakb.start <= peaka.start - margin <= peakb.end or peakb.start <= peaka.end - margin <= peakb.end:
            return True
        if peakb.start <= peaka.start + margin <= peakb.end or peakb.start <= peaka.end + margin <= peakb.end:
            return True
        if peaka.start <= peakb.start - margin <= peaka.end or peaka.start <= peakb.end - margin <= peaka.end:
            return True
        if peaka.start <= peakb.start + margin <= peaka.end or peaka.start <= peakb.end + margin <= peaka.end:
            return True
    return False


def stichPeaks(peaks, margin=1):
    """
    Stich close peaks based on postion array. Peaks are all in the same chromosome
    """
    cov = set()
    for i, peak in enumerate(peaks):
        cov.update(range(peak.start, peak.end + 1))
    cov = list(cov)
    cov.sort()
    npeaks = []
    i = 0
    while i < len(cov) - 1:
        for j in range(i + 1, len(cov)):
            if cov[j] - cov[j - 1] > margin:
                break
            else:
                continue
        #print(i,j,len(cov[i:j]),cov[j-1],cov[j]) #used for debug
        #peak = peaks[0].chrom
        peak = Peak()
        peak.chrom = peaks[0].chrom
        peak.start = cov[i]
        peak.end = cov[j - 1]
        peak.length = cov[j - 1] - cov[i] + 1
        npeaks.append(peak)
        i = j  #update search start
    return npeaks


#loops related operations
def checkAnchorOverlap(xa, xb, ya, yb):
    """
    check the overlap of a region for the same chromosome
    """
    if (ya <= xa <= yb) or (ya <= xb <= yb) or (ya <= xa <= xb <= yb):
        return True
    if (xa <= ya <= xb) or (xa <= yb <= xb) or (xa <= ya <= yb <= xb):
        return True
    return False


def checkLoopOverlap(loopa, loopb):
    """
    check the overlap of two loops
    """
    if loopa.chromX != loopb.chromX or loopa.chromY != loopb.chromY:
        return False
    if checkAnchorOverlap(loopa.x_start, loopa.x_end, loopb.x_start,
                          loopb.x_end) and checkAnchorOverlap(
                              loopa.y_start, loopa.y_end, loopb.y_start,
                              loopb.y_end):
        return True
    return False


def combineLoops(loops, loops_2):
    """
    Combine loops result.
    """
    for key in loops_2.keys():
        if key not in loops:
            loops[key] = loops_2[key]
        else:
            ds = set()
            for loop in loops[key]:
                r = [
                    loop.chromX, loop.x_start, loop.x_end, loop.chromY,
                    loop.y_start, loop.y_end
                ]
                ds.add(tuple(r))
            for loop in loops_2[key]:
                r = [
                    loop.chromX, loop.x_start, loop.x_end, loop.chromY,
                    loop.y_start, loop.y_end
                ]
                if tuple(r) not in ds:
                    loops[key].append(loop)
    return loops
