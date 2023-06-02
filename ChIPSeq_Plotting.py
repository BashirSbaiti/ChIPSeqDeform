import os

import pandas as pd
import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt
import matplotlib as mpl
import regex as re
import dill

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
LEN_CHR = 4699742

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

### Plot number of peaks in expanding complementary sequence (fig 4d,f)

if PLOT_EXPANDING_SEED:

    def countPeaksInSpan(peakSer, span):
        """ :param peakSer: series of peaks where index = location, value = intensity
            :param span: array like of (inclusive, exclusive) specifying span to count peaks in"""
        return len([val for ind, val in peakSer.items() if span[0] <= ind < span[1]])


    # binding site for all gRNAs = gRNA seq + NGG
    AGGAAseq = "GGGACCTAAGATTTGAGGAA" + ".GG"
    ACCCAseq = "TCGAACACACTCTCTACCCA" + ".GG"
    LacZseq = "TCGTTTTACAACGTCGTGAC" + ".GG"

    bindSites = [LacZseq, AGGAAseq, ACCCAseq]
    names = ["lacZ", "AGGAA", "ACCCA"]

    ### import referance genome in the same way as main.py
    RefGenomeFname = "Rostain Supp Info/MG1655 Referance Genome/GCA_000005845.2_ASM584v2_genomic.fna"
    record = SeqIO.parse(RefGenomeFname, "fasta")
    refObj = list(record)[0]
    refGenome = (refObj.seq, refObj.seq.reverse_complement())

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
        plt.close(fig)
        sampleInd += 1

#### Now we investigate deformability stuff

## Plot binding activity over the genome as a rolling average

window = 100000
# list(arr[slice]) + list(arr) copies the last window-1 elements again, so they are not wiped out by rolling function
allPeaksDfExpand = pd.concat([pd.Series(list(arr[-(window-1):]) + list(arr), name=name) for arr, name in zip(allArrs, names)], axis=1)
raAllPeaksDfExpand = allPeaksDfExpand.rolling(window, axis=0).mean()
rstdAllPeaksDfExpand = allPeaksDfExpand.rolling(window, axis=0).std()
raAllPeaksDf = raAllPeaksDfExpand.iloc[window-1:,:].reset_index(drop=True)
rstdAllPeaksDf = rstdAllPeaksDfExpand.iloc[window-1:,:].reset_index(drop=True)

fig, axs = plt.subplots(len(names), 1, sharex="all", sharey="all", figsize=(8, 7))
colors = ["red", "green", "blue"]

for colName, ax, color in zip(raAllPeaksDf.columns, axs, colors):
    rollAvg = raAllPeaksDf.loc[:, colName]
    rollStd = rstdAllPeaksDf.loc[:, colName]

    ax.plot(rollAvg, label=colName, color=color)
    ax.legend()
    ax.set_ylabel("Intensity")
    ax.fill_between(np.arange(0, len(rollAvg), 1), rollAvg - 0.0*rollStd, rollAvg + 0.0*rollStd, color=color, alpha=0.5)

fig.suptitle(f"All Binding Events Across the Genome for all gRNA, Averaged w/Surrounding {window}bp")
plt.xticks(range(0, len(allArrs[0]), len(allArrs[0]) // 24), rotation=45)
axs[-1].set_xlabel("Location on MG1655 Genome")
plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig(f"figures/AveragedBindingAcrossGenomeAll")
