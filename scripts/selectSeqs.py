import numpy as np
from Bio import SeqIO
import pandas as pd

def selectSequences(fromFile, toFile, selCrit=lambda i,rec:True):
    fromCount = 0
    toCount = 0
    with open(toFile, 'w') as w:
        for i, rec in enumerate(SeqIO.parse(fromFile, "fasta")):
            fromCount += 1
            if selCrit(i, rec): 
                toCount += 1
                #print(f">{rec.description}\n{rec.seq[65:-75]}", file=w)
                print(f">{rec.description}\n{''.join(['N']*65)}{str(rec.seq[65:-75])}{''.join(['N']*75)}", file=w)
    print (f'Wrote {toCount}/{fromCount} seqs to {toFile}, replacing N-/C-termini with N')
fastafile = "../results/aligned_Gisaid0525_b1_b2_presel.fasta"
fastafile = "aligned_2020-06-25_b1-3QCNR_countryMax100_k25.fasta"
fastafileOut = "aligned_2020-06-25_b1-3QCNR_countryMax100_k25O.fasta"

prun=['Turkey/6224-Ankara1034/2020', 'Canada/ON_MU-S53/2020', 'SouthAfrica/R02606/2020',
      'BosniaandHerzegovina/03_Tuzla/2020'] 

#selCrit = lambda i,rec: rec.description.split('/')[0] in ['UAE', 'UnitedArabEmirates', 'Wuhan-Hu-1'] 
selCrit = lambda i,rec: not rec.id in prun
selectSequences(fastafile, fastafileOut, selCrit=selCrit)
