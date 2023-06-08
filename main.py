import matplotlib.pyplot as plt
import numpy as np
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import GEOparse
import urllib.request
import gzip
import shutil
import regex as re
import pandas as pd
import dill

GEOparse.logger.set_verbosity("ERROR")

IMPORTFILES = False

if IMPORTFILES:

    ############# Data Importing ###################

    Entrez.email = "bs407@duke.edu"


    def searchDB(db, searchTerm):
        handle = Entrez.esearch(db=db, term=searchTerm, retmax=200)
        record = Entrez.read(handle)
        ids = record['IdList']
        return ids


    #### ChIP-seq data #####

    # Search for ChIP-Seq data on GEO
    ids = searchDB("gds", "GSE217004")
    print(f"Found {len(ids)} hits: {ids}.")

    # There are 11 hits:
    # header file, raw data all together, 3 * triplicates gRNAs (1 targeting, 2 nontargeting (toxic)) processed (9 total)

    # remove header and raw data
    ids.remove("200217004")
    ids.remove("100021222")

    # Summarize processed ChIP Seq Data
    accessionNums = list()
    for id in ids:
        handle = Entrez.esummary(db="gds", id=id, report="full")
        record = Entrez.read(handle)
        accessionNums.append(record[0]["Accession"])

    # Retrieve data, save to machine
    for acc in accessionNums:
        gse = GEOparse.get_GEO(geo=acc, destdir='GEOFiles')
        suppFTPAdd = gse.metadata['supplementary_file_1'][0]
        suppFTPAdd = suppFTPAdd.replace("fftp://", "ftp://")  # idk what fftp is but im gonna pretend it doesnt exist
        name = suppFTPAdd[suppFTPAdd.rfind("/") + 1:]
        urllib.request.urlretrieve(suppFTPAdd, f"GEOFiles/{name}")

        # unzip downloaded file
        with gzip.open(f"GEOFiles/{name}", "rb") as fIn:
            with open(f"GEOFiles/{name}".replace(".gz", ""), "wb") as fOut:
                shutil.copyfileobj(fIn, fOut)

    ##### Referance Genome #####

    # Search for Genome data
    ids = searchDB("assembly", "ASM584v2")
    print(f"Found {len(ids)} hits: {ids}.")

    # There is only one hit
    id = ids[0]

    # find data file
    handle = Entrez.esummary(db="assembly", id=id, report="full")
    record = Entrez.read(handle)
    FTPUrl = record['DocumentSummarySet']["DocumentSummary"][0]['FtpPath_GenBank']
    FileName = "GCA_000005845.2_ASM584v2_genomic.fna.gz"
    FTPUrl += f"/{FileName}"

    # save file
    name = FTPUrl[FTPUrl.rfind("/") + 1:]
    urllib.request.urlretrieve(FTPUrl, f"Rostain Supp Info/MG1655 Referance Genome/{FileName}")

    # unzip downloaded file
    with gzip.open(f"Rostain Supp Info/MG1655 Referance Genome/{FileName}", "rb") as fIn:
        with open(f"Rostain Supp Info/MG1655 Referance Genome/{FileName}".replace(".gz", ""), "wb") as fOut:
            shutil.copyfileobj(fIn, fOut)

    ############# Main Body ###################

RefGenomeFname = "Rostain Supp Info/MG1655 Referance Genome/GCA_000005845.2_ASM584v2_genomic.fna"

record = SeqIO.parse(RefGenomeFname, "fasta")  # iterable object with metadata about seq and finally seq object with seq
refObj = list(record)[0]  # there is only one seq in this file (e coli only have one chromosome)
# ALL SEQUENCES GO 5 TO 3

refGenome = (refObj.seq, refObj.seq.reverse_complement())


# 5 to 3 forward and 5 to 3 reverse strands


class GenomicTgtSite:

    def __init__(self, refGenome, PAMLoc, contextLen=500):
        """ :param refGenome: tuple representing reference genome (5-3 fwd strand, 5-3 reverse strand)
            :param PAMLoc: location of PAM site in question tuple of (strand num, 5' start incl, 3' stop excl)
            :param contextLen: number of bp to pull on either side of PAM+tgt site for local context"""
        self.GENOME = refGenome
        self.PAMLoc = PAMLoc
        self.gRNA = self.getgRNA()
        self.context5, self.context3 = self.getContext(contextLen) if self.getContext(contextLen) is not None else (None, None)
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
        else: # e. coli genome is ciruclar
            if self.PAMLoc[1] - 20 >= 0:
                context5Fw = fwStrand[context5Start:] + fwStrand[0:self.PAMLoc[1] - 20]
            else:
                context5Fw = fwStrand[context5Start:self.PAMLoc[1]-20]
        context5Rv = context5Fw.reverse_complement()

        if len(context5Fw) > length:
            raise Exception(f"5' is at {len(strand)}")

        context3End = self.PAMLoc[2] + length
        if context3End <= len(fwStrand) - 1:
            context3Fw = fwStrand[self.PAMLoc[2]:context3End]
        else: # e. coli genome is ciruclar
            overflow = context3End - len(fwStrand)
            context3Fw = fwStrand[self.PAMLoc[2]:] + fwStrand[0:overflow]
        context3Rv = context3Fw.reverse_complement()

        if len(context3Fw) > length:
            raise Exception(f"5' is at {len(strand)}")

        return (context5Fw, context5Rv), (context3Fw, context3Rv)

    def getDeform(self):
        """ :returns tuple of (avg deform L, avg deform R) over context"""

        deformRefDF = pd.read_csv("DeformabilityParameters_AlbertoPerez.txt", delimiter="\t", index_col=0).T

        strand = self.context5[0]
        if len(strand) >= 2:
            deformList = list()  # list of pd series containing deform params for every dinuc in the seq
            for dinucStart in range(0, len(strand)-1):
                dinucStop = dinucStart+2
                dinuc = strand[dinucStart:dinucStop]
                dinuc = dinuc.reverse_complement() if str(dinuc) not in deformRefDF.columns else dinuc
                dinucDeform = deformRefDF.loc[:, str(dinuc)]
                deformList.append(dinucDeform)
            context5DeformDf = pd.concat(deformList, axis=1)
            context5DeformAvg = context5DeformDf.mean(axis=1)
        else:
            context5DeformAvg = None

        # now do it again for 3'
        strand = self.context3[0]
        if len(strand) >= 2:
            deformList = list()
            for dinucStart in range(0, len(strand)-1):
                dinucStop = dinucStart+2
                dinuc = strand[dinucStart:dinucStop]
                dinuc = dinuc.reverse_complement() if str(dinuc) not in deformRefDF.columns else dinuc
                dinucDeform = deformRefDF.loc[:, str(dinuc)]
                deformList.append(dinucDeform)
            context3DeformDf = pd.concat(deformList, axis=1)
            context3DeformAvg = context3DeformDf.mean(axis=1)
        else:
            context3DeformAvg = None

        return context5DeformAvg, context3DeformAvg


allPAMIndxs = list()

for strandNum, strand in enumerate(refGenome):
    for val in re.finditer(r".GG", str(strand), overlapped=True):
        allPAMIndxs.append(val.span()[0] if strandNum == 0 else len(strand) - 1 - val.span()[0])

plt.figure()
plt.hist(allPAMIndxs, bins=1000)

allPAMScores = np.zeros(len(refGenome[0]))
for i in allPAMIndxs:
    allPAMScores[i] += 1

scoreWindow = 2
allScoresDfExpand = pd.Series(list(allPAMScores[-(scoreWindow - 1):]) + list(allPAMScores), name="Markers at all PAMs")
raAllScoresExpand = allScoresDfExpand.rolling(scoreWindow, axis=0).mean()
rstdAllScoresExpand = allScoresDfExpand.rolling(scoreWindow, axis=0).std()
raAllScores = raAllScoresExpand.iloc[scoreWindow - 1:].reset_index(drop=True)
rstdAllScores = rstdAllScoresExpand.iloc[scoreWindow - 1:].reset_index(drop=True)

plt.figure()
plt.plot(raAllScores)
plt.title(scoreWindow)
plt.show()


# TODO: calculating Tm could be a good control
