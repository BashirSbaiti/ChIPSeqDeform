import os

import pandas as pd
import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt
import matplotlib as mpl
import regex as re
from Bio.Seq import Seq
import dill
import time

plt.tight_layout()

LEN_CHR = 4641652
PLOT_GENOME_DEFORM = False


def getStrandDeform(strand):
    deformRefDF = pd.read_csv("DeformabilityParameters_AlbertoPerez.txt", delimiter="\t", index_col=0).T

    if len(strand) >= 2:
        deformList = list()  # list of pd series containing deform params for every dinuc in the seq
        for dinucStart in range(0, len(strand) - 1):
            dinucStop = dinucStart + 2
            dinuc = strand[dinucStart:dinucStop]
            dinuc = dinuc.reverse_complement() if str(dinuc) not in deformRefDF.columns else dinuc
            dinucDeform = deformRefDF.loc[:, str(dinuc)]
            deformList.append(dinucDeform)
        deformDf = pd.concat(deformList, axis=1)
        deformDfAvg = deformDf.mean(axis=1)
    else:
        deformDfAvg = None

    return deformDfAvg


class GenomicTgtSite:

    def __init__(self, refGenome, PAMLoc, contextLen=500):
        """ :param refGenome: tuple representing reference genome (5-3 fwd strand, 5-3 reverse strand)
            :param PAMLoc: location of PAM site in question tuple of (strand num, 5' start incl, 3' stop excl)
            :param contextLen: number of bp to pull on either side of PAM+tgt site for local context"""
        self.GENOME = refGenome
        self.PAMLoc = PAMLoc
        self.gRNA = self.getgRNA()
        self.context5, self.context3 = self.getContext(contextLen) if self.getContext(contextLen) is not None else (
            None, None)
        self.deform5, self.deform3 = self.getDeform()

    def __repr__(self):
        PRINT_DEFORM = True
        VERBOSE_SEQ = False
        if PRINT_DEFORM:
            str = f"Target site with PAMLoc = {self.PAMLoc}; gRNA = {self.gRNA}; Context5' = {self.context5}" \
                  f" (len = {len(self.context5[0])} bp, deform = {self.deform5}); Context3' = {self.context3} (len = {len(self.context3[0])} bp, deform = {self.deform3})"
        else:
            str = f"Target site with PAMLoc = {self.PAMLoc}; gRNA = {self.gRNA}; Context5' = {self.context5}" \
                  f" (len = {len(self.context5[0])} bp); Context3' = {self.context3} (len = {len(self.context3[0])})"
            if VERBOSE_SEQ:
                str = f"Target site with PAMLoc = {self.PAMLoc}; gRNA = {self.gRNA}; Context5' = {self.context5[0]}" \
                      f" (len = {len(self.context5[0])} bp); Context3' = {self.context3[0]} (len = {len(self.context3[0])} bp)"
        return str

    def getgRNA(self):
        """ finds target site sequence given PAM location (read 5-3) on the genome"""
        strand = self.GENOME[self.PAMLoc[0]]
        tsStart = self.PAMLoc[1] - 20
        if tsStart >= 0:
            ret = Seq(str(strand[self.PAMLoc[1] - 20:self.PAMLoc[1]]).replace("T", "U"))
        else:
            ret = Seq(str(strand[tsStart:] + strand[0:self.PAMLoc[1]]).replace("T", "U"))
        return ret

    def getContext(self, length):
        """ Find (length) bp fw and rev (fw is defined as PAM containing) on each side of target site"""

        fwStrand = self.GENOME[self.PAMLoc[0]]

        context5Start = self.PAMLoc[1] - 20 - length
        if context5Start >= 0:
            context5Fw = fwStrand[context5Start:self.PAMLoc[1] - 20]
        else:  # e. coli genome is ciruclar
            if self.PAMLoc[1] - 20 >= 0:
                context5Fw = fwStrand[context5Start:] + fwStrand[0:self.PAMLoc[1] - 20]
            else:
                context5Fw = fwStrand[context5Start:self.PAMLoc[1] - 20]
        context5Rv = context5Fw.reverse_complement()

        if len(context5Fw) > length:
            raise Exception(f"5' is at {len(context5Fw)}")

        context3End = self.PAMLoc[2] + length
        if context3End <= len(fwStrand) - 1:
            context3Fw = fwStrand[self.PAMLoc[2]:context3End]
        else:  # e. coli genome is ciruclar
            overflow = context3End - len(fwStrand)
            context3Fw = fwStrand[self.PAMLoc[2]:] + fwStrand[0:overflow]
        context3Rv = context3Fw.reverse_complement()

        if len(context3Fw) > length:
            raise Exception(f"5' is at {len(context3Fw)}")

        return (context5Fw, context5Rv), (context3Fw, context3Rv)

    def getDeform(self):
        """ :returns tuple of (avg deform L, avg deform R) over context"""
        context5DeformAvg, context3DeformAvg = getStrandDeform(self.context5[0]), getStrandDeform(self.context3[0])

        return context5DeformAvg, context3DeformAvg


### import referance genome in the same way as main.py
RefGenomeFname = "Rostain Supp Info/MG1655 Referance Genome/GCA_000005845.2_ASM584v2_genomic.fna"
record = SeqIO.parse(RefGenomeFname, "fasta")
refObj = list(record)[0]
refGenome = (refObj.seq, refObj.seq.reverse_complement())

with open("allPeaksDfReps.pkl", "rb") as dill_file:
    allPeaksDfReps = dill.load(dill_file)

# Looking at deform across the genome
# We will use nonoverlapping chunks of chunkSize bp
deformAxs = ["Twist", "Tilt", "Roll", "Shift", "Slide", "Rise"]

if PLOT_GENOME_DEFORM:
    chunkSize = 1000
    allSers = list()
    for chunkStart in range(0, LEN_CHR - (chunkSize) + 1, chunkSize):
        chunkStop = chunkStart + chunkSize
        chunk = refGenome[0][chunkStart:chunkStop]
        chunkDeform = pd.Series(getStrandDeform(chunk), name=chunkStart, dtype=float)
        allSers.append(chunkDeform)
    allChunksDeformDf = pd.concat(allSers, axis=1)

    for ax in deformAxs:
        deformOfAx = allChunksDeformDf.loc[ax, :]
        plt.figure()
        plt.title(f"Distribution of Deform in {ax}, Genome Chunks of {chunkSize}")
        plt.boxplot(deformOfAx.array)
        # plt.yticks(np.linspace(plt.ylim()[0], plt.ylim()[1], 20))
        # plt.hist(deformOfAx.array, bins=1000)
        # plt.savefig(f"figures/deform/{ax}DistAcrossGenomeChunk{chunkSize}")
        plt.savefig(f"figures/deform/{ax}DistAcrossGenomeChunk{chunkSize}Box")
        plt.close()

# Redoing expanding seed cals

AGGAAseq = "GGGACCTAAGATTTGAGGAA" + ".GG"
ACCCAseq = "TCGAACACACTCTCTACCCA" + ".GG"
LacZseq = "TCGTTTTACAACGTCGTGAC" + ".GG"

bindSites = [LacZseq, AGGAAseq, ACCCAseq]
names = ["lacZ", "AGGAA", "ACCCA"]

deformAxis = "Rise"
medDeformLimits = (7.78, 7.80)
# inclusive on both sides

for name, bindSite in zip(names, bindSites):
    peakSers = allPeaksDfReps.loc[:, name]
    numReps = len(peakSers.columns)
    seqs = [bindSite[start:] for start in range(len(bindSite) - 1 - 9, len(bindSite) - 2)]
    fracSers = list()  # index = deform (LMH), value = fracBound/total for each rep, one per overlap seq
    sitesSers = list()  # index = deform (LMH), value = total for each rep, one per overlap seq
    avgIntensitySers = list()  # index = deform (LMH), value = average intensity of all bind events (across diff reps), one per overlap seq
    stdIntensitySers = list()
    # series of differing complementarity to gRNA+PAM
    count = 0
    for seq in seqs:
        peaksPerDeform = {"Low": np.zeros(numReps), "Medium": np.zeros(numReps), "High": np.zeros(numReps)}  # key = deform category (LMH), value is list of num peaks for each rep
        sitesPerDeform = {"Low": np.zeros(numReps), "Medium": np.zeros(numReps), "High": np.zeros(numReps)}  # key = LMH, val = list of totalSites for each rep
        intensitiesPerDeform = {"Low": [], "Medium": [], "High": []}  # key = LMH, val = list of intensities across all reps
        #count += 1
        #if count == 3:
        #    break
        # only using fw strand of referance strand because I assume that if they find binding to the reverse strand they map it to the fw
        for find in re.finditer(seq, str(refGenome[0]), overlapped=False):
            ts = GenomicTgtSite(refGenome, (0, find.span()[1] - 3, find.span()[1]), contextLen=500)
            deform = float(ts.deform5[deformAxis] + ts.deform5[deformAxis]) / 2.0
            if medDeformLimits[0] <= deform <= medDeformLimits[1]:
                deformCat = "Medium"
            elif deform > medDeformLimits[1]:
                deformCat = "High"
            else:
                deformCat = "Low"
            for repNum in range(numReps): # for each replicate of each gRNA
                peakSer = peakSers.iloc[:, repNum]
                sitePeaks = 0
                siteIntensity = 0
                for pos in range(find.span()[0], find.span()[1]):  # search for a peak at all positions
                    if peakSer.loc[pos] > 0:
                        sitePeaks = 1
                        siteIntensity = peakSer.loc[pos]
                        break  # do not double count peaks

                if sitePeaks == 0:  # if no peaks are found, the intensity at this site is zero
                    siteIntensity = 0

                peaksPerDeform[deformCat][repNum] += sitePeaks
                intensitiesPerDeform[deformCat].append(siteIntensity)
                sitesPerDeform[deformCat][repNum] += 1

        sitesSer = pd.Series(sitesPerDeform, name=seq)
        sitesSers.append(sitesSer)

        avgIntensityPerDeform = {}
        stdIntensityPerDeform = {}
        for deformCat in intensitiesPerDeform.keys():
            avgIntensityPerDeform[deformCat] = (np.sum(intensitiesPerDeform[deformCat]) / np.sum(sitesPerDeform[deformCat])) if \
                np.sum(sitesPerDeform[deformCat]) != 0 else 0.0
            stdIntensityPerDeform[deformCat] = np.std(intensitiesPerDeform[deformCat]) if len(
                intensitiesPerDeform[deformCat]) > 0 else 0.0
        avgIntensitySer = pd.Series(avgIntensityPerDeform, name=seq)
        stdIntensitySer = pd.Series(stdIntensityPerDeform, name=seq)
        avgIntensitySers.append(avgIntensitySer)
        stdIntensitySers.append(stdIntensitySer)

        fracsPerDeform = {}
        for deformCat in sitesPerDeform.keys():
            fracsPerDeform[deformCat] = [(peaksPerDeform[deformCat][repNum] / sitesPerDeform[deformCat][repNum]) if sitesPerDeform[
                                                                                                       deformCat][repNum] != 0 else 0.0 for repNum in range(numReps)]
        fracSer = pd.Series(fracsPerDeform, name=seq)
        fracSers.append(fracSer)

    def elementWiseAvg(element):
        return np.average(element)
    def elementWiseStd(element):
        return np.std(element)

    # plot deform seperated line chart with fracs
    deformSepFracDf = pd.concat(fracSers, axis=1) # each element is a list of len 3 for 3 reps
    deformSepFracDfAvg = deformSepFracDf.apply(np.vectorize(elementWiseAvg))
    deformSepFracDfStd = deformSepFracDf.apply(np.vectorize(elementWiseStd))
    fig = plt.figure(figsize=(8, 6))
    deformSepFracDfAvg.T.plot(yerr=deformSepFracDfStd.T, linestyle='dashed', marker='o', ms=3)
    plt.title(f"Fraction Bound Sites for {name}")
    plt.ylabel("Fraction of Sites With a Peak")
    plt.xticks(rotation=45)
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(f"figures/deform/DeformExpandingSeed{name}Reps")
    plt.close()

    # plot bar chart of sites per deform
    deformSepStiesDf = pd.concat(sitesSers, axis=1)
    deformSepStiesDfAvg = deformSepStiesDf.apply(np.vectorize(elementWiseAvg))
    deformSepStiesDfStd = deformSepStiesDf.apply(np.vectorize(elementWiseStd))
    fig = plt.figure()
    deformSepStiesDfAvg.T.plot(kind="bar", logy=True, yerr=deformSepStiesDfStd.T)
    plt.title(f"Number of Sites By Deform for {name}")
    plt.ylabel("Number of Sites")
    plt.xticks(rotation=45)
    plt.subplots_adjust(bottom=0.25)
    plt.savefig(f"figures/deform/DeformExpandingSeedCountSites{name}Reps")
    plt.close()

    # plot bar chart of avg intensity per deform
    deformSepIntenAvgDf = pd.concat(avgIntensitySers, axis=1)
    deformSepIntenStdDf = pd.concat(stdIntensitySers, axis=1)
    fig = plt.figure()
    deformSepIntenAvgDf.T.plot(kind="bar", yerr=deformSepIntenStdDf.T, logy=True)
    plt.title(f"Average Intensity for Peaks By Deform for {name}")
    plt.ylabel("Average Peak Intensity")
    plt.xticks(rotation=45)
    plt.subplots_adjust(bottom=0.25)
    plt.savefig(f"figures/deform/DeformExpandingSeedIntensity{name}Reps")
    plt.close()
