import os
import time

import pandas as pd
import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt
import regex as re
from Bio.Seq import Seq

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
        self.deform = self.getDeform()

    def __repr__(self):
        PRINT_DEFORM = True
        VERBOSE_SEQ = False
        if PRINT_DEFORM:
            str = f"Target site with PAMLoc = {self.PAMLoc}; gRNA = {self.gRNA}; Context5' = {self.context5}" \
                  f" (len = {len(self.context5[0])} bp; Context3' = {self.context3} (len = {len(self.context3[0])} bp; deform = {self.deform})"
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
        """ :returns series of average deform over context for all axes"""
        context5DeformAvg, context3DeformAvg = getStrandDeform(self.context5[0]), getStrandDeform(self.context3[0])

        return (context5DeformAvg + context3DeformAvg) / 2


### import referance genome in the same way as file_importing.py
RefGenomeFname = "Rostain Supp Info/MG1655 Referance Genome/GCA_000005845.2_ASM584v2_genomic.fna"
record = SeqIO.parse(RefGenomeFname, "fasta")
refObj = list(record)[0]
refGenome = (refObj.seq, refObj.seq.reverse_complement())

if __name__ == "__main__":

    # with open("allPeaksDfReps.pkl", "rb") as dill_file:
    #    allPeaksDfReps = dill.load(dill_file)

    allPeaksDfReps = pd.read_csv("allPeaksDfReps.csv", sep="\t", encoding="utf-8", index_col=0)
    # for some reason pandas renames duplicate columns names when reading csv's, below code will reverse this renaming
    for colname in allPeaksDfReps.columns:
        allPeaksDfReps.rename(columns={colname: colname if colname.find(".") == -1 else colname[0:colname.find(".")]},
                              inplace=True)

    # Looking at deform across the genome
    # We will use nonoverlapping chunks of chunkSize bp
    deformAxs = ["Twist", "Tilt", "Roll", "Shift", "Slide", "Rise"]

    if not os.path.isfile("allChunksDeformDf.csv"):  # does file exist yet
        chunkSize = 1000
        allSers = list()
        for chunkStart in range(0, LEN_CHR - (chunkSize) + 1, chunkSize):
            chunkStop = chunkStart + chunkSize
            chunk = refGenome[0][chunkStart:chunkStop]
            chunkDeform = pd.Series(getStrandDeform(chunk), name=chunkStart, dtype=float)
            allSers.append(chunkDeform)
        allChunksDeformDf = pd.concat(allSers, axis=1)
        allChunksDeformDf.to_csv("allChunksDeformDf.csv", sep="\t", encoding="utf-8")
    else:
        allChunksDeformDf = pd.read_csv("allChunksDeformDf.csv", sep="\t", encoding="utf-8", index_col=0)
    if PLOT_GENOME_DEFORM:
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

    # Redoing expanding seed calcs with deform

    AGGAAseq = "GGGACCTAAGATTTGAGGAA" + ".GG"
    ACCCAseq = "TCGAACACACTCTCTACCCA" + ".GG"
    LacZseq = "TCGTTTTACAACGTCGTGAC" + ".GG"

    bindSites = [LacZseq, AGGAAseq, ACCCAseq]
    names = ["lacZ", "AGGAA", "ACCCA"]

    deformAxis = "Rise"
    fracsHL = [.18]  # top what frac and bottom what frac are considered high and low deform respectively
    medDeformLimitss = list()  # list of tuples defining hard numeric cutoffs for med

    # calculate hard numeric cutoffs based on fracsHL
    deformOfAx = allChunksDeformDf.loc[deformAxis, :].sort_values()
    for fracHL in fracsHL:
        lowCutoff = deformOfAx.iloc[int(fracHL * len(deformOfAx))]
        highCutoff = deformOfAx.iloc[int((1 - fracHL) * len(deformOfAx))]
        medDeformLimitss.append((lowCutoff, highCutoff))

    # load precalculated deforms (all 500bp context)
    # TODO: label files with 500bp for support with other context lens
    start = time.perf_counter()
    precalcSeqToDf = {}  # key = precalculated seq, value = df of deform at all occurances for all axes
    for fileName in os.listdir():
        if fileName.find("allGenomic") != -1 and fileName.find("Deform.csv") != -1:
            seq = fileName[fileName.find("Genomic")+len("Genomic"):fileName.find("Deform")]
            seqRe = seq.replace("N", ".")

            df = pd.read_csv(fileName, delimiter="\t", index_col=0)

            precalcSeqToDf[seq] = df

    end = time.perf_counter()
    print(f"finished precalc file importing in {end-start} sec")
    for medDeformLimits, fracHL in zip(medDeformLimitss, fracsHL):
        for name, bindSite in zip(names, bindSites):
            peakSers = allPeaksDfReps.loc[:, name]
            numReps = len(peakSers.columns)
            seqs = [bindSite[start:] for start in
                    range(len(bindSite) - 1 - 9, len(bindSite) - 2)]  # series of differing complementarity to gRNA+PAM
            fracSers = list()  # index = deform (LMH), value = fracBound/total for each rep, one per overlap seq
            sitesSers = list()  # index = deform (LMH), value = total for each rep, one per overlap seq
            intensitySers = list()  # index = deform (LMH), value = intensity all peaks across reps and in same rep, one per overlap seq

            foundPAMEnds = list()  # keep track of end of PAM of every seq we find so we dont double count one ...ngg in multiple categories

            count = 0
            for seq in seqs:
                peaksPerDeform = {"Low": np.zeros(numReps), "Medium": np.zeros(numReps), "High": np.zeros(
                    numReps)}  # key = deform category (LMH), value is list of num peaks for each rep
                sitesPerDeform = {"Low": np.zeros(numReps), "Medium": np.zeros(numReps),
                                  "High": np.zeros(numReps)}  # key = LMH, val = list of totalSites for each rep
                intensitiesPerDeform = {"Low": [], "Medium": [],
                                        "High": []}  # key = LMH, val = list of intensities across all reps
                # count += 1
                # if count == 4:
                #     break
                startTime = time.perf_counter()
                # only using fw strand of referance strand because I assume that if they find binding to the reverse strand they map it to the fw
                for find in re.finditer(seq, str(refGenome[0]), overlapped=False):
                    if find.span()[1] not in foundPAMEnds:
                        if seq != ".GG": # we search genome for ngg's last, so we only need to make sure identified ngg isnt part of a larger match
                            foundPAMEnds.append(find.span()[1])

                        # calculate deform
                        if seq in precalcSeqToDf.keys():  # if it was already calculated, no need to do it again!
                            deformDf = precalcSeqToDf[seq]
                            deform = deformDf.loc[deformAxis, find.span()[0]]
                        else:
                            ts = GenomicTgtSite(refGenome, (0, find.span()[1] - 3, find.span()[1]), contextLen=500)
                            deform = ts.deform[deformAxis]
                        if medDeformLimits[0] <= deform <= medDeformLimits[1]:
                            deformCat = "Medium"
                        elif deform > medDeformLimits[1]:
                            deformCat = "High"
                        else:
                            deformCat = "Low"
                        for repNum in range(numReps):  # for each replicate of each gRNA
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
                            if siteIntensity != 0:
                                intensitiesPerDeform[deformCat].append(siteIntensity)
                            sitesPerDeform[deformCat][repNum] += 1
                endTime = time.perf_counter()
                print(f"Finished {seq} in {endTime-startTime} sec")
                intensitySer = pd.Series(intensitiesPerDeform, name=seq)
                intensitySers.append(intensitySer)

                sitesSer = pd.Series(sitesPerDeform, name=seq)
                sitesSers.append(sitesSer)

                fracsPerDeform = {}
                for deformCat in sitesPerDeform.keys():
                    fracsPerDeform[deformCat] = [
                        (peaksPerDeform[deformCat][repNum] / sitesPerDeform[deformCat][repNum]) if sitesPerDeform[
                                                                                                       deformCat][
                                                                                                       repNum] != 0 else 0.0
                        for repNum in range(numReps)]
                fracSer = pd.Series(fracsPerDeform, name=seq)
                fracSers.append(fracSer)

            deformSepFracDf = pd.concat(fracSers, axis=1)  # each element is a list of len 3 for 3 reps

            deformSepStiesDf = pd.concat(sitesSers, axis=1)

            deformSepIntenDf = pd.concat(intensitySers,
                                         axis=1)  # each element is a list variable len as each val corresponda to 1 peak
            # because this would make a very big pkl file, just save sum, std, avg for plotting instead of whole list
            deformSepIntenSumDf = deformSepIntenDf.apply(np.vectorize(np.sum))
            deformSepIntenAvgDf = deformSepIntenDf.apply(np.vectorize(np.average))
            deformSepIntenStdDf = deformSepIntenDf.apply(np.vectorize(np.std))

            dfsToSave = {"deformSepFracDf": deformSepFracDf,
                         "deformSepStiesDf": deformSepStiesDf,
                         "deformSepIntenSumDf": deformSepIntenSumDf,
                         "deformSepIntenAvgDf": deformSepIntenAvgDf,
                         "deformSepIntenStdDf": deformSepIntenStdDf
                         }
            destinationDir = f"DeformPlottingDfs/{deformAxis}/Cutoff={medDeformLimits}={fracHL * 100}%/{name}"
            if not os.path.exists(destinationDir):
                os.makedirs(destinationDir)
            for varName in dfsToSave.keys():
                # with open(f"{destinationDir}/{varName}.pkl", "wb") as dill_file:
                #    dill.dump(dfsToSave[varName], dill_file)
                dfsToSave[varName].to_csv(f"{destinationDir}/{varName}.csv", sep="\t", encoding="utf-8")
