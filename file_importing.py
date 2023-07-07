import matplotlib.pyplot as plt
import numpy as np
from Bio import Entrez
from Bio import SeqIO
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


allPAMIndxs = list()

for strandNum, strand in enumerate(refGenome):
    for val in re.finditer(r".GG", str(strand), overlapped=True):
        pass # found a pam site


# TODO: calculating Tm could be a good control
