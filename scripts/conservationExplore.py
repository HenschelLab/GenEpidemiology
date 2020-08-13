from collections import Counter
from Bio import SeqIO
import pdb
from scipy.stats import entropy 


def codons(rec, S1, pos1, pos2 ):
    try:
        sq = rec[S1 - 1 + pos1*3 : S1 - 1 + pos2*3]
        aa = sq.translate()
    except:
        return None
        #pdb.set_trace()
    return str(aa.seq)

def u(rec):
    return rec.description.split('/')[0] in ['UAE', 'UnitedArabEmirates']

cdss = [[266, 13468],[13468, 21555],[21563, 25384],[25393, 26220],[26245, 26472],[26523, 27191]]

S1,S2 = 21563,25384 ## feature of S protein according to reference.gb for Wuhan-Hu sequence (in config dir)
furinPos1, furinPos2 = 680, 686
mutS1, mutS2 = 613, 614
fasta = "../results/aligned_Gisaid0525_b1_b2_presel.fasta"

#codonStats = Counter([codons(rec, S1, furinPos1, furinPos2) for rec in SeqIO.parse(fasta, "fasta") if u(rec)])
codonStats = Counter([codons(rec, S1, furinPos1, furinPos2) for rec in SeqIO.parse(fasta, "fasta")])
#codonStats = Counter([codons(rec, S1, mutS1, mutS2) for rec in SeqIO.parse(fasta, "fasta") if u(rec)])
