import pandas as pd
import numpy as np
from Bio import SeqIO
import regex as re
from Bio.Seq import Seq
from ChIPSeq_Deform import GenomicTgtSite, refGenome

seq = "NGG"
seqRe = seq.replace("N", ".")  # regex standard formatting

LEN_CHR = len(refGenome[0])

allTsSers = list()

for strandNum, strand in enumerate(refGenome):
    for val in re.finditer(seqRe, str(strand), overlapped=True):
        dict = {}
        ts = GenomicTgtSite(refGenome, (strandNum, val.span()[0], val.span()[1]), contextLen=500)
        dict["PAM Location"] = val.span()[0] if strandNum == 0 else LEN_CHR - 1 - val.span()[0]
        for deformAx in ts.deform.index:
            dict[deformAx] = ts.deform[deformAx]
        tsSer = pd.Series(dict)
        allTsSers.append(tsSer)

allTsDf = pd.concat(allTsSers, axis=1)
allTsDf.to_csv(f"allGenomic{seq}Deform.csv", sep="\t", encoding="utf-8", header=False)
