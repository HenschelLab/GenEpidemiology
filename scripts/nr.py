## Redundancy reduction by country
# run ../scripts/nr.py aligned_2020-06-25_b1-3QC.fasta aligned_2020-06-25_b1-3QC_NR.fasta

import hashlib
from Bio import SeqIO
import sys, os
from collections import defaultdict


def nr(infile, outfile, checkExist=True):
    totalCounter = 0
    nrCount = 0
    nrDict = defaultdict(set)

    if checkExist and os.path.exists(outfile):
        print("Warning: {outfile} exists, aborting")
        sys.exit(1)

    with open(outfile, 'w') as w:
        for rec in SeqIO.parse(infile, 'fasta'):
            totalCounter += 1
            country = rec.id.split('/')[0]
            fingerprint = hashlib.sha224(str(rec.seq[65:-75]).encode('ascii')).hexdigest()
            if country in ['UAE','UnitedArabEmirates'] or not fingerprint in nrDict[country]:
                nrCount += 1
                print(f'>{rec.id}\n{rec.seq}', file=w)
                nrDict[country].add(fingerprint)

    print(f"Reduction: {totalCounter} => {nrCount} ({100*nrCount/totalCounter:.2f}%)")

if __name__ == "__main__": 
    infile, outfile = sys.argv[-2:]
    nr(infile, outfile)
