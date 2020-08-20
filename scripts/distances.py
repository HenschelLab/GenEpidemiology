# run preprocess.py ../results/aligned_iDownLocal2.fasta

from Bio import SeqIO
import numpy as np
import sys

code = {'A': '0001', 'C':'0010', 'G':'0100', 'T':'1000', 'N':'0000'}  
tt = "ACTGN".maketrans(code) 
mask = {'A': '0000', 'C':'0000', 'G':'0000', 'T': '0000', 'N':'1111'}
ttm = "ACTGN".maketrans(mask) 

class BinRep:
    def __init__(self, rec):
        self.fromUAE = rec.description.split('/')[0] in ['UAE', 'UnitedArabEmirates']
        self.dna = ''.join([nt if nt in set('ACGTN') else 'N' for nt in rec.seq.upper()])
        self.binRep = self.dna.translate(tt)
        self.mask = self.dna.translate(ttm)
        self.a = int(self.binRep, 2)
        self.m = int(self.mask, 2)

    def bindiff(self, other): ## roughly 10x speedup in comparison to previous seq. distance count
        return bin((self.a^other.a) & ~(self.m|other.m)).count('1')//2

if __name__ == "__main__":
    seqFile = sys.argv[-1]
    all2 = [BinRep(rec) for rec in SeqIO.parse(seqFile, "fasta")] 
    uae2 = [seq for seq in all2 if seq.fromUAE]
    d = np.zeros((len(all2),len(uae2)))

    ## pairwise distances
    for i, seq in enumerate(all2):
        for j, uaeseq in enumerate(uae2):
            d[i,j] = uaeseq.bindiff(seq)
    np.save('../results/distances_b3a', d)
