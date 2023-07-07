import os

import pandas as pd
import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt
import matplotlib as mpl
import regex as re
import dill
import time

plt.tight_layout()

PLOT_EXPANDING_SEED = False
ON_CLUSTER = False

if not ON_CLUSTER:
    import seaborn as sns

if ON_CLUSTER:
    mpl.use("Agg")
else:
    mpl.use("Qt5Agg")


class Peak:
    def __init__(self, location, intensity):
        """:param location: tuple of [start, stop)
           :param intensity: intensity of the peak (float)"""
        self.location = location
        self.intensity = intensity

    def __repr__(self):
        return f"ChIP-Seq peak at {self.location} w/intensity {self.intensity}"

    def getLoc(self):
        return self.location

    def getIntensity(self):
        return self.intensity


CHROMOSOME = "NC_000913"
LEN_CHR = 4641652

allPeakSers = list()  # list of series, each series contains all peaks for one ChIP seq file

for file in os.listdir("GEOFiles"):
    if file.find(".bed") != -1 and file.find(".gz") == -1:
        name = file[file.find("ChIP") + 5:file.find("summits") - 1]
        peaks = list()
        with open(f"GEOFiles/{file}", 'r') as f:
            for line in f:
                line = line.replace("\n", "")
                lst = line.split("\t")
                chr = lst[0]
                loc = (int(lst[1]), int(lst[2]))
                intensity = float(lst[4])
                if chr == CHROMOSOME:
                    ipeak = Peak(loc, intensity)
                    peaks.append(ipeak)

        ser = pd.Series(peaks, name=name, dtype=object)
        allPeakSers.append(ser)

allPeakArrLabeled = list()
for ser in allPeakSers:
    peakArr = np.zeros(LEN_CHR)
    for el in ser:
        peakArr[el.getLoc()[0]] = el.getIntensity()
    peakArrLabeled = pd.Series(peakArr, name=ser.name[:ser.name.find("-")])
    allPeakArrLabeled.append(peakArrLabeled)
peaksArrDf = pd.concat(allPeakArrLabeled, axis=1)

#with open("allPeaksDfReps.pkl", "wb") as dill_file:
#    dill.dump(peaksArrDf, dill_file)

peaksArrDf.to_csv("allPeaksDfReps.csv", sep="\t", encoding="utf-8")

    ## Verify that all loc distances are 1, so we can just use the start loc as peak pos
for ser in allPeakSers:
    for el in ser:
        loc = el.getLoc()
        if (loc[1] - loc[0]) != 1:
            raise Exception(f"loc difference is {loc[1] - loc[0]} not 1")

### Plot number of binding events per gRNA

names = list()
binds = list()
for ser in allPeakSers:
    names.append(ser.name[0:ser.name.find("-")])
    binds.append(len(ser))

serNames = pd.Series(names, name="name")
serBinds = pd.Series(binds, name="bind events")
dfEventsPergRNA = pd.concat([serNames, serBinds], axis=1)

if not ON_CLUSTER:
    sns.barplot(data=dfEventsPergRNA, x="name", y="bind events", errorbar=('ci', 90))
    plt.savefig(f"figures/BindingEventsPergRNA")
    plt.close()

# Turn series into arrays where index = start loc, intensity = val at index
# load all these arrays into a list and add replicates

allArrs = list()

for ser in allPeakSers:
    if ser.name.find("-1") != -1:
        arr = np.zeros(LEN_CHR, dtype=float)
    else:
        arr = allArrs[-1]

    for peak in ser:
        arr[peak.getLoc()[0]] += peak.getIntensity()

    if ser.name.find("-1") != -1:
        allArrs.append(arr)

allArrs = [arr / 3 for arr in allArrs]
# len = 3 [lacZ, AGGAA, ACCCA]
names = ["lacZ", "AGGAA", "ACCCA"]

for grnaNum in range(len(allArrs)):
    for i, el in enumerate(allArrs[grnaNum]):
        if el != 0:
            start = i - 10
            stop = i + 10
            for j, surr in enumerate(allArrs[grnaNum][start:stop]):
                if surr != 0 and i != start+j:
                    allArrs[grnaNum][i] += surr
                    allArrs[grnaNum][start+j] = 0

allArrsNonzero = list()
# nonzero version of allArrs, each list element still corresponds to one gRNA but instead represented by a pd series
# where index = peak location, data = intensity

for i, arr in enumerate(allArrs):
    x = [ind for ind in range(len(arr)) if arr[ind] != 0]
    y = [arr[ind] for ind in range(len(arr)) if arr[ind] != 0]
    allArrsNonzero.append(pd.Series(data=y, index=x, name=names[i]))

### plot distribution of peak intensities for all gRNAs
for ser in allArrsNonzero:
    plt.figure()
    plt.hist([val for val in ser], bins=160)
    plt.ylabel("Count")
    plt.xlabel("Intensity")
    plt.title(f"Distribution of Peak Intensities for {ser.name}")
    plt.savefig(f"figures/DistIntensities{ser.name}")
    plt.close()

### plot binding events across a window for all gRNAs

# window = (3723500, 3727500) # for fig 4a
window = (3882500, 3886500)  # for fig 4b
# window = (0, LEN_CHR)
PlotLargest10All = False
# (inclusive, exclusive)

for i, ser in enumerate(allArrsNonzero):
    x = [ind for ind in ser.index if window[0] <= ind < window[1]]  # all nonzero peak positions in window
    y = [val for ind, val in ser.items() if window[0] <= int(ind) < window[1]]  # all nonzero intensities in window
    if PlotLargest10All:
        largest10 = ser.nlargest(10)
        print(f"Largest 10 binding events for {ser.name}: {largest10}")
    fig = plt.figure(figsize=(8, 6))
    plt.title(f"All Binding Events Across the Genome for gRNA {names[i]}")
    plt.xlabel("Location on MG1655 Genome")
    plt.ylabel("Intensity")
    plt.plot(x, y)
    plt.xticks(range(window[0], window[1], (window[1] - window[0]) // 24), rotation=45)
    if PlotLargest10All:
        plt.vlines(largest10.index, 0, largest10.values, colors="r")
    plt.tight_layout()
    plt.savefig(f"figures/AllBindingAcrossGenome{names[i]}")
    # plt.show()
    plt.close()

### Plot number of peaks in expanding complementary sequence (fig 4d,f)

# binding site for all gRNAs = gRNA seq + NGG
AGGAAseq = "GGGACCTAAGATTTGAGGAA" + ".GG"
ACCCAseq = "TCGAACACACTCTCTACCCA" + ".GG"
LacZseq = "TCGTTTTACAACGTCGTGAC" + ".GG"

bindSites = [LacZseq, AGGAAseq, ACCCAseq]
names = ["lacZ", "AGGAA", "ACCCA"]

### import referance genome in the same way as file_importing.py
RefGenomeFname = "Rostain Supp Info/MG1655 Referance Genome/GCA_000005845.2_ASM584v2_genomic.fna"
record = SeqIO.parse(RefGenomeFname, "fasta")
refObj = list(record)[0]
refGenome = (refObj.seq, refObj.seq.reverse_complement())

if PLOT_EXPANDING_SEED:

    def countPeaksInSpan(peakSer, span):
        """ :param peakSer: series of peaks where index = location, value = intensity
            :param span: array like of (inclusive, exclusive) specifying span to count peaks in"""
        return 1 if len([val for ind, val in peakSer.items() if span[0] <= ind < span[1]]) > 0 else 0


    sampleInd = 0
    for name, bindSite in zip(names, bindSites):
        seqs = [bindSite[start:] for start in range(len(bindSite) - 1 - 9, len(bindSite) - 2)]
        fracs = list()
        # series of differing complementarity to gRNA+PAM
        for seq in seqs:
            totalPeaks = 0
            totalSites = 0
            # only using fw strand of referance strand because I assume that if they find binding to the reverse strand they map it to the fw
            for val in re.finditer(seq, str(refGenome[0]), overlapped=True):
                peaks = countPeaksInSpan(allArrsNonzero[sampleInd], val.span())
                totalPeaks += peaks
                totalSites += 1
            frac = totalPeaks / totalSites
            fracs.append(frac)
            print(f"{seq} in {name}: {totalPeaks} peaks, {totalSites} sites, {frac} peaks per site")
        fig = plt.figure(figsize=(8, 6))
        plt.plot(seqs, fracs, linestyle="dashed")
        plt.ylabel("Fraction of Sites With a Peak")
        plt.xticks(rotation=45)
        plt.subplots_adjust(bottom=0.15)
        plt.savefig(f"figures/ExpandingSeed{name}")
        plt.close()
        sampleInd += 1

## Plot binding activity over the genome as a rolling average

bindWindow = 10000
allPeaksDf = pd.concat([pd.Series(arr, name=name) for arr, name in zip(allArrs, names)],
                       axis=1)  # this will not be useful unless we use window = 1 = no RA (used in later calculations)
# list(arr[slice]) + list(arr) copies the last bindWindow-1 elements again, so they are not wiped out by rolling function
allPeaksDfExpand = pd.concat(
    [pd.Series(list(arr[-(bindWindow - 1):]) + list(arr), name=name) for arr, name in zip(allArrs, names)], axis=1)
raAllPeaksDfExpand = allPeaksDfExpand.rolling(bindWindow, axis=0).mean()
rstdAllPeaksDfExpand = allPeaksDfExpand.rolling(bindWindow, axis=0).std()
raAllPeaksDf = raAllPeaksDfExpand.iloc[bindWindow - 1:, :].reset_index(drop=True)
rstdAllPeaksDf = rstdAllPeaksDfExpand.iloc[bindWindow - 1:, :].reset_index(drop=True)

fig, axs = plt.subplots(len(names), 1, sharex="all", sharey="all", figsize=(8, 7))
colors = ["red", "green", "blue"]

for colName, ax, color in zip(raAllPeaksDf.columns, axs, colors):
    rollAvg = raAllPeaksDf.loc[:, colName]
    rollStd = rstdAllPeaksDf.loc[:, colName]

    ax.plot(rollAvg, label=colName, color=color)
    ax.legend()
    ax.set_ylabel("Intensity")
    ax.fill_between(np.arange(0, len(rollAvg), 1), rollAvg - 0.0 * rollStd, rollAvg + 0.0 * rollStd, color=color,
                    alpha=0.5)

fig.suptitle(f"All Binding Events Across the Genome for all gRNA, Averaged w/Surrounding {bindWindow}bp")
plt.xticks(range(0, len(allArrs[0]), len(allArrs[0]) // 24), rotation=45)
axs[-1].set_xlabel("Location on MG1655 Genome")
plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig(f"figures/AveragedBindingAcrossGenomeAll")
plt.close()

print(f"{time.time()}: Completed All Before Score Calcs!")

## Impletment plotting of gRNA overlap and see how that varies across genome


allOverlapDicts = list()  # key = overlaping sequence, value = length of overlap
for ts, name in zip(bindSites, names):
    overlapDict = {
        ".GG": 3}  # we define PAM is 3 to differentiate it from no seq overlap with no PAM site (this way overlap
                   # len+ngg = int), we will later subtract 3 everywhere to get gRNA seq overlap
    for start in range(len(ts) - 4, -1, -1):
        overlapDict[ts[start:]] = len(ts[start:])
    allOverlapDicts.append(overlapDict)

allOverlapsArr = list()
for overlapDict, name in zip(allOverlapDicts, names):
    overlapsArr = np.zeros(LEN_CHR)
    for strandNum, strand in enumerate(refGenome):
        strand = str(strand)
        for val in re.finditer(r".GG", strand, overlapped=True):
            start, stop = val.span()
            exceededMaxMatch = False
            while not exceededMaxMatch:
                start -= 1
                if start < 0:
                    raise Exception("Invalid start index when calculating gRNA overlap. Implement circularization.")
                newSeq = strand[start:stop]
                generalNewSeq = newSeq[:-3] + ".GG"
                if generalNewSeq not in overlapDict.keys():
                    exceededMaxMatch = True
            start += 1
            match = strand[start:stop]
            generalMatch = match[:-3] + ".GG"
            overlaps = overlapDict[generalMatch]
            if strandNum == 0:
                newStart = start
                newStop = stop
            elif strandNum == 1:
                newStart = len(overlapsArr) - stop
                newStop = len(overlapsArr) - start
            overlapsArr[newStart:newStop] = [max(overlaps - 2, val) for val in overlapsArr[newStart:newStop]]
    allOverlapsArr.append(overlapsArr)

print(f"{time.time()}: Completed gRNA Overlap Calculations!")
# first, take rolling avg
for overlapWindow in [10000]:
    allOverlapsDf = pd.concat([pd.Series(arr, name=name) for arr, name in zip(allOverlapsArr, names)],
                              axis=1)  # make no rollavg version for use later
    if overlapWindow == 1:  # if overlapWindow == 1, don't do rolling average
        raAllOverlapsDf = allOverlapsDf
        rstdAllOverlapsDf = 0 * raAllOverlapsDf
    else:
        # list(arr[slice]) + list(arr) copies the last overlapWindow-1 elements again, so they are not wiped out by rolling function
        allOverlapsDfExpand = pd.concat(
            [pd.Series(list(arr[-(overlapWindow - 1):]) + list(arr), name=name) for arr, name in zip(allOverlapsArr, names)],
            axis=1)
        raAllOverlapsDfExpand = allOverlapsDfExpand.rolling(overlapWindow, axis=0).mean()
        rstdAllOverlapsDfExpand = allOverlapsDfExpand.rolling(overlapWindow, axis=0).std()
        raAllOverlapsDf = raAllOverlapsDfExpand.iloc[overlapWindow - 1:, :].reset_index(drop=True)
        rstdAllOverlapsDf = rstdAllOverlapsDfExpand.iloc[overlapWindow - 1:, :].reset_index(drop=True)

    print(f"{time.time()}: Completed Overlap Rolling Avg!")
    # then, plot
    fig, axs = plt.subplots(len(names), 1, sharex="all", sharey="all", figsize=(8, 7))
    colors = ["red", "green", "blue"]

    for colName, ax, color in zip(raAllOverlapsDf.columns, axs, colors):
        rollAvg = raAllOverlapsDf.loc[:, colName]
        rollStd = rstdAllOverlapsDf.loc[:, colName]

        ax.scatter(np.arange(0, len(rollAvg), 1), rollAvg, s=3, alpha=0.5, label=colName, color=color)
        ax.legend()
        ax.set_ylabel("gRNA Overlap")
        # ax.fill_between(np.arange(0, len(rollAvg), 1), rollAvg - 0.0 * rollStd, rollAvg + 0.0 * rollStd, color=color,
        # alpha=0.5)

    fig.suptitle(f"gRNA Overlap across the genome for all gRNA, rolling window = {overlapWindow}bp")
    axs[-1].set_xlabel("Location on MG1655 Genome")
    plt.ylabel("gRNA Overlap")
    plt.xticks(range(0, LEN_CHR, (LEN_CHR) // 24), rotation=45)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(f"figures/AllOverlapsComparisonRA{overlapWindow}")
    plt.close()

## Plot overlaps vs binding activity side by side to see if there is a relation

for name in names:
    rollAvgOverlap = raAllOverlapsDf.loc[:, name]
    rollStdOverlap = rstdAllOverlapsDf.loc[:, name]

    rollAvgBind = raAllPeaksDf.loc[:, name]
    rollStdBind = rstdAllPeaksDf.loc[:, name]

    fig, axs = plt.subplots(2, 1, sharex="all", figsize=(8, 5))
    fig.suptitle(f"gRNA Overlap vs Bind Activity for {name}, Rolling Avg {overlapWindow}bp")
    if overlapWindow != bindWindow:
        raise Exception(f"gRNA Overlap rolling avg window {overlapWindow} is not equal to bind window {bindWindow}.")

    axs[0].plot(rollAvgOverlap, label=f"{name} gRNA Overlap")
    axs[0].legend()
    axs[0].set_ylabel("gRNA Overlap")

    axs[1].plot(rollAvgBind, label=f"{name} Binding")
    axs[1].legend()
    axs[1].set_ylabel("Binding Intensity")
    axs[1].set_xlabel("Location on MG1655 Genome")

    plt.xticks(range(0, LEN_CHR, LEN_CHR // 24), rotation=45)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(f"figures/OverlapvsBinding{name}SideBySide")
    plt.close()

# Same as above, but plot one against the other
USE_ROLLAVG = False

# TODO: we dont currently take into account 0bp (pam only binding events), may be worth looking into

for name in names:
    if USE_ROLLAVG:
        overlaps = raAllOverlapsDf.loc[:, name]
        stdOverlaps = rstdAllOverlapsDf.loc[:, name]

        bind = raAllPeaksDf.loc[:, name]
        stdBind = rstdAllPeaksDf.loc[:, name]
    else:
        overlaps = allOverlapsDf.loc[:, name]
        stdOverlaps = 0 * overlaps

        bind = allPeaksDf.loc[:, name]
        stdBind = 0 * bind

    fig = plt.figure()
    if USE_ROLLAVG:
        plt.title(f"Binding Intensity vs gRNA Overlap for {name}, Rolling Avg {overlapWindow}bp")
    else:
        plt.title(f"Binding Intensity vs gRNA Overlap for {name}")

    if overlapWindow != bindWindow and USE_ROLLAVG:
        raise Exception(f"gRNA Overlap rolling avg window {overlapWindow} is not equal to bind window {bindWindow}.")

    # filter out zero intensity peaks
    scoreFiltered = overlaps[bind != 0]
    bindFiltered = bind[bind != 0]

    plt.scatter(x=scoreFiltered, y=bindFiltered, s=3)
    plt.xlabel("gRNA Overlap")
    plt.ylabel("ChIP-Seq Binding Intensity")
    plt.xlim(left=-0.3)
    plt.xticks(np.arange(0, max(scoreFiltered) + 1, 1))

    plt.tight_layout()
    if USE_ROLLAVG:
        plt.savefig(
            f"figures/ScorevsBinding{name}OneAgainstOtherRA{bindWindow}")
    else:
        plt.savefig(
            f"figures/OverlapvsBind{name}.svg")
    plt.close()

# Save allPeaks, allBinds info to use in other files

#with open("allPeaksDf.pkl", "wb") as dill_file:
#    dill.dump(allPeaksDf, dill_file)

allPeaksDf.to_csv("allPeaksDf.csv", sep="\t", encoding="utf-8")

#with open("allOverlapsDf.pkl", "wb") as dill_file:
#    dill.dump(allOverlapsDf, dill_file)

allOverlapsDf.to_csv("allOverlapsDf.csv", sep="\t", encoding="utf-8")
