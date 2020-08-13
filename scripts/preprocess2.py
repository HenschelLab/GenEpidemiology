import numpy as np
from Bio import SeqIO
import pandas as pd
import hashlib
from collections import defaultdict


chineseCities = [ 'Anhui',  'Beijing', 'Chongqing', 'Fuzhou', 'Foshan', 'Fujian', 'Fuyang', 'Ganzhou','Guangdong', 'Guangzhou', 'Hangzhou', 'Harbin', 'Hefei', 'Henan', 'HongKong', 
                  'Jian', 'Jiangsu', 'Jiangxi', 'Jingzhou', 'Jiujiang',  'Lishui',  'MN908947_S',  'NanChang', 'Nanchang',
                  'Pingxiang',  'Shandong', 'Shaoxing', 'Shanghai', 'Shangrao', 'Shenzhen', 'Sichuan',   'Tianmen', 'Wuhan', 'Wuhan-Hu-1', 'Xinyu', 'Yichun', 'Yingtan', 'Yunnan', 'Zhejiang']

                  #'canine', 'cat', 'env', 'hCoV-19_S', 'mink', 'tiger'

#colors = [cc.cm.CET_D1A(c) for c in np.linspace(0,1.05, maxScore)] ##WTF goes wrong here: NoneType error

def nr(infile, outfile):
    nrDict = defaultdict(set)
    repDict = {}
    keepers = []

    for rec in SeqIO.parse(infile, 'fasta'):
        country = dfi.loc[rec.id]['country']
        # hash from reliable 
        fingerprint = hashlib.sha224(str(rec.seq[65:-75]).encode('ascii')).hexdigest()
        if country in ['UAE','UnitedArabEmirates'] or not fingerprint in nrDict[country]:
            nrCount += 1
            nrDict[country].add(fingerprint)
            repDict[fingerprint] = rec.id
            keepers.append(rec.id)
            #dfi.loc[rec.id]['repres'] = repDict.get(fingerprint, rec.id)
    return keepers
            

def getcountry(loc):
    if loc in chineseCities: return 'China'
    elif loc == 'UnitedArabEmirates': return 'UAE'
    elif loc == 'ITALY': return 'Italy'
    elif loc in 'England Scotland Wales NorthernIreland'.split(): return 'UK'
    else: return loc

def selectSequences(fromFile, toFile, selCrit=lambda i,rec:True):
    counter = 0
    with open(toFile, 'w') as w:
        for i, rec in enumerate(SeqIO.parse(fromFile, "fasta")):
            if selCrit(i, rec): 
                if rec.id == 'Wuhan/Hu-1/2019': rec.id = 'MN908947'
                counter += 1
                #print(f">{rec.description}\n{rec.seq}", file=w)
                ## replaceing tail and end with N's (so to unify possibly low qual regions)
                print(f">{rec.description}\n{''.join(['N']*65)}{str(rec.seq[65:-75])}{''.join(['N']*75)}", file=w)
    print(f"Wrote {counter} sequences to {toFile}")

def seqStats(rec):
    return (rec.description, getcountry(rec.description.split('/')[0]),
            len(rec.seq), rec.seq.count('-'), rec.seq.count('N'), rec.seq[65:-73].count('N'))

k=25
countryMax = 50

def downsample(g, country, alreadySel, notyetSel):
    g0 = alreadySel[alreadySel.country==country]
    freeSpots = countryMax - len(g0)
    if freeSpots <= 0: return set()
    eligible = notyetSel[notyetSel.country==country]
    spots = min(freeSpots, len(eligible))
    result = set(eligible.sample(spots).index)  
#    if len(result) < min(len(eligible), countryMax):
#        import pdb; pdb.set_trace()
    return result 

#fastafile = "../data/sequences_2020-06-25_b1-3QC.fasta"
fastafile = "../results/aligned_2020-06-25_b1-3QC.fasta"
fastafileNR = "../results/aligned_2020-06-25_b1-3QC_NR.fasta"
toFile1 = f"../results/aligned_2020-06-25_b1-3QCNR_k{k}.fasta"
toFile2 = f"../results/aligned_2020-06-25_b1-3QCNR_countryMax{countryMax}_k{k}.fasta"
toFile3 = f"../results/aligned_2020-06-25_b1-3QCNR_countryMax{countryMax}_k{k}_1k.fasta"
#fastafile = "../results/aligned_2020-06-25_b1b2.fasta"

uae = [rec.name for rec in SeqIO.parse(fastafile, "fasta") if rec.description.split('/')[0] in ['UAE', 'UnitedArabEmirates']]
df = pd.DataFrame([seqStats(rec) for rec in SeqIO.parse(fastafile, "fasta")]) 
df.columns = 'strain country length gaps N Nclipped'.split() 
d = np.load('../results/distances_b3a.npy')
ddf = pd.DataFrame(d, index=df.strain, columns=uae)
dfi = df.set_index('strain')
dfi['repres'] = None


if True:
    ## NR
    keepers = nr(fastafile, fastafileNR)
    ## Downsampling
    ddfM = pd.concat([dfi, ddf], axis=1, sort=False)
    ddfMQ = ddfM[ddfM.Nclipped<30].loc[keepers]
    
    ## choose k-NN for each emirati sample, 
    s = set([])
    columns = ddfMQ.columns.values[6:]
    for col in columns:
        s = s.union(ddfMQ[ddfMQ[col]<10].sort_values(by=col).index.values[:k])

    ## downsample each country (friggin UK all over the place)
    alreadySel = ddfMQ.loc[s]
    notyetSel  = ddfMQ.loc[set(ddfMQ.index)-s]
    subsample = [downsample(g, c, alreadySel, notyetSel) for c, g in ddfMQ.groupby('country')]  #ddfMQ.loc[s] for limiting to UAE-close seqs
    s1 = set(s).union(*subsample) 
    
#    nohup augur tree --alignment aligned_2020-06-25_b1-3QCNR_countryMax300.fasta --output tree_raw_2020-06-25_b1-3QCNR_countryMax300.nwk --nthreads 60

    #s1 = s1.union(ddfMQ[ddfMQ.country.isin(selectedLocations)].sample(50).index.values)  
    s1 = s1.union(uae) ## shouldn't be necessary
    s1.add('Wuhan/Hu-1/2019') ## they changed the naming 
    #s_1k = set(ddfMQ.loc[s1-s].sample(1000).index).union(s) 
    #s_1k = s_1k.union(uae) ## shouldn't be necessary
    #s_1k.add('Wuhan/Hu-1/2019') ## they changed the naming 

    ## UAE-close sequences
    selCrit1 = lambda i,rec: rec.name in s
    ## UAE-close + random (with country limit)
    selCrit2 = lambda i,rec: rec.name in s1
    ## UAE-close + 1000 random (with country limit)    
    selCrit3 = lambda i,rec: rec.name in s_1k

#    selectSequences(fromFile=fastafile, toFile=toFile1, selCrit=selCrit1)
    selectSequences(fromFile=fastafile, toFile=toFile2, selCrit=selCrit2)
#    selectSequences(fromFile=fastafile, toFile=toFile3, selCrit=selCrit3)


"""
[('India', 200),
 ('Australia', 200),
 ('Denmark', 200),
 ('Iceland', 200),
 ('UK', 200),
 ('USA', 200),
 ('Spain', 200),
 ('Sweden', 200),
 ('China', 167),
 ('Russia', 154),
 ('Belgium', 142),
 ('Switzerland', 124),
 ('NewZealand', 124),
 ('Canada', 104),
 ('SouthAfrica', 104),
 ('Germany', 103),
 

In [30]: ddf[ddf.index.str.startswith('UAE')]                                                                                                                                                                                               
Out[30]: 
         UAE/16B  UAE/12B  UAE/13  UAE/16  UAE/14B  UAE/2  UAE/31B
UAE/16B      0.0      8.0     9.0    10.0     11.0   12.0     17.0
UAE/12B      8.0      0.0     1.0     4.0      3.0    6.0     11.0
UAE/13       9.0      1.0     0.0     5.0      4.0    7.0     12.0
UAE/16      10.0      4.0     5.0     0.0      5.0    6.0     13.0
UAE/14B     11.0      3.0     4.0     5.0      0.0    9.0     14.0
UAE/2       12.0      6.0     7.0     6.0      9.0    0.0     15.0
UAE/31B     14.0      8.0     9.0    10.0     11.0   12.0      0.0

In [31]: ddf[ddf.index.str.startswith('UnitedArab')]                                                                                                                                                                                        
Out[31]: 
                               UAE/16B  UAE/12B  UAE/13  UAE/16  UAE/14B  UAE/2  UAE/31B
UnitedArabEmirates/L0484/2020     12.0      6.0     7.0     8.0      9.0   10.0      9.0
UnitedArabEmirates/L1758/2020     12.0      6.0     7.0     8.0      9.0   10.0      9.0
UnitedArabEmirates/L2185/2020     12.0      6.0     7.0     8.0      9.0   10.0      9.0
UnitedArabEmirates/L0/2020        15.0      9.0    10.0    11.0     12.0   13.0     12.0
UnitedArabEmirates/L5621/2020     11.0      5.0     6.0     7.0      8.0    9.0      8.0
UnitedArabEmirates/L7356/2020     12.0      6.0     7.0     8.0      9.0   10.0      9.0
UnitedArabEmirates/L9440/2020     14.0      8.0     9.0    10.0     11.0   12.0     11.0
UnitedArabEmirates/L1076/2020     15.0      9.0    10.0    11.0     12.0   13.0     12.0
UnitedArabEmirates/L9768/2020     10.0      4.0     5.0     6.0      7.0    8.0     13.0
UnitedArabEmirates/L0881/2020      7.0      1.0     2.0     3.0      4.0    5.0     10.0
UnitedArabEmirates/L3779/2020      8.0      2.0     3.0     4.0      5.0    6.0     11.0
UnitedArabEmirates/L0879/2020     10.0      4.0     5.0     6.0      7.0    8.0     13.0
UnitedArabEmirates/L5630/2020     18.0     12.0    13.0    14.0     15.0   16.0     15.0
UnitedArabEmirates/L4280/2020     15.0      9.0    10.0    11.0     12.0   13.0     12.0
UnitedArabEmirates/L1014/2020     15.0      9.0    10.0    11.0     12.0   13.0     12.0
UnitedArabEmirates/L4184/2020     16.0     10.0    11.0    12.0     13.0   14.0     13.0
UnitedArabEmirates/L9766/2020     16.0     10.0    11.0    12.0     13.0   14.0     13.0
UnitedArabEmirates/L6599/2020     20.0     14.0    15.0    16.0     17.0   18.0     17.0
UnitedArabEmirates/L2409/2020     16.0     12.0    13.0    14.0     15.0   16.0     15.0
UnitedArabEmirates/L6627/2020     17.0     13.0    14.0    15.0     16.0   17.0     16.0
UnitedArabEmirates/L4682/2020     19.0     13.0    14.0    15.0     16.0   17.0     16.0
UnitedArabEmirates/L0184/2020     21.0     17.0    18.0    19.0     20.0   21.0     20.0
UnitedArabEmirates/L0904/2020     22.0     18.0    19.0    20.0     21.0   22.0     21.0
UnitedArabEmirates/L0231/2020     23.0     17.0    18.0    19.0     20.0   21.0     20.0
UnitedArabEmirates/L068/2020      21.0     17.0    18.0    19.0     20.0   21.0     20.0


In [61]: mind=ddfMQ.loc[s]                                                                                                              

In [62]: mind['minD'] = mind[columns].min(axis=1)                                                                                       

In [63]: mind                                                                                                                           
Out[63]: 
                               country  length  gaps    N  ...  UAE/H12_S14  UnitedArabEmirates/L6841/2020  UAE/H29_S21  minD
strain                                                     ...                                                               
Canada/Qc-L00241344/2020        Canada   29903     0   46  ...          5.0                            7.0          6.0   3.0
Spain/Granada-COV002889/2020     Spain   29903     0   87  ...         11.0                           11.0         12.0   2.0
England/20132054704/2020       England   29903     0   70  ...          8.0                            8.0          9.0   3.0
Wuhan/HB-WH1-140/2020            China   29903     0    0  ...          8.0                            8.0          9.0   3.0
England/201040079/2020         England   29903     0   58  ...         10.0                           10.0         11.0   1.0
...                                ...     ...   ...  ...  ...          ...                            ...          ...   ...
Scotland/EDB024/2020          Scotland   29903     0  129  ...         10.0                           10.0         11.0   1.0
Belgium/ITM_C067/2020          Belgium   29903     0  134  ...         14.0                           14.0         15.0   5.0
England/CAMB-73D5C/2020        England   29903     0  124  ...         13.0                           13.0         14.0   1.0
Sichuan/SC-WCH-082/2020          China   29903     0    0  ...          7.0                            7.0          8.0   2.0
Wuhan/HB-WHCM-104/2020           China   29903     0    0  ...          8.0                           10.0          9.0   4.0

[952 rows x 63 columns]

In [64]: mind[mind.minD<3]                                                                                                              
Out[64]: 
                               country  length  gaps    N  ...  UAE/H12_S14  UnitedArabEmirates/L6841/2020  UAE/H29_S21  minD
strain                                                     ...                                                               
Spain/Granada-COV002889/2020     Spain   29903     0   87  ...         11.0                           11.0         12.0   2.0
England/201040079/2020         England   29903     0   58  ...         10.0                           10.0         11.0   1.0
Scotland/EDB021/2020          Scotland   29903     0  123  ...         10.0                           10.0         11.0   1.0
England/20142025204/2020       England   29903     0   47  ...         13.0                           13.0         14.0   1.0
England/20132081104/2020       England   29903     0   52  ...         10.0                           10.0         11.0   1.0
...                                ...     ...   ...  ...  ...          ...                            ...          ...   ...
Iceland/516/2020               Iceland   29903     0   45  ...          6.0                            6.0          7.0   1.0
Germany/BAV-MVP0039/2020       Germany   29903     0  128  ...         10.0                           10.0         11.0   1.0
Scotland/EDB024/2020          Scotland   29903     0  129  ...         10.0                           10.0         11.0   1.0
England/CAMB-73D5C/2020        England   29903     0  124  ...         13.0                           13.0         14.0   1.0
Sichuan/SC-WCH-082/2020          China   29903     0    0  ...          7.0                            7.0          8.0   2.0

[507 rows x 63 columns]

In [66]: Counter(mind[mind.minD<3].country)                                                                                             
Out[66]: 
Counter({'Spain': 34,
         'England': 156,
         'Scotland': 18,
         'Switzerland': 13,
         'China': 32,
         'Australia': 16,
         'India': 10,
         'Russia': 16,
         'UAE': 57,
         'USA': 46,
         'Belgium': 5,
         'Oman': 9,
         'Denmark': 10,
         'Canada': 4,
         'Sweden': 3,
         'Iceland': 17,
         'SouthAfrica': 9,
         'Vietnam': 3,
         'Wales': 19,
         'Germany': 7,
         'Indonesia': 2,
         'Jordan': 1,
         'NewZealand': 10,
         'Turkey': 1,
         'Italy': 2,
         'Qatar': 2,
         'Chile': 2,
         'Colombia': 1,
         'France': 1,
         'MN908947': 1})

In [68]: Counter(mind[mind.minD<3].country).most_common()                                                                               
Out[68]: 
[('England', 156),
 ('UAE', 57),
 ('USA', 46),
 ('Spain', 34),
 ('China', 32),
 ('Wales', 19),
 ('Scotland', 18),
 ('Iceland', 17),
 ('Australia', 16),
 ('Russia', 16),
 ('Switzerland', 13),
 ('India', 10),
 ('Denmark', 10),
 ('NewZealand', 10),
 ('Oman', 9),
 ('SouthAfrica', 9),
 ('Germany', 7),
 ('Belgium', 5),
 ('Canada', 4),
 ('Sweden', 3),
 ('Vietnam', 3),
 ('Indonesia', 2),
 ('Italy', 2),
 ('Qatar', 2),
 ('Chile', 2),
 ('Jordan', 1),
 ('Turkey', 1),
 ('Colombia', 1),
 ('France', 1),
 ('MN908947', 1)]

In [69]: Counter(mind[mind.minD<2].country).most_common()                                                                               
Out[69]: 
[('England', 123),
 ('UAE', 57),
 ('USA', 26),
 ('Wales', 15),
 ('Iceland', 12),
 ('Russia', 9),
 ('Spain', 9),
 ('SouthAfrica', 9),
 ('Scotland', 8),
 ('Switzerland', 8),
 ('Australia', 7),
 ('India', 7),
 ('Germany', 6),
 ('NewZealand', 6),
 ('Oman', 4),
 ('Denmark', 4),
 ('Canada', 4),
 ('Vietnam', 3),
 ('Sweden', 2),
 ('Italy', 2),
 ('Qatar', 2),
 ('China', 2),
 ('Belgium', 1),
 ('France', 1),
 ('MN908947', 1)]

In [71]: Counter(mind[mind.minD<1].country).most_common()                                                                               
Out[71]: [('UAE', 57), ('Canada', 1)]

In [72]: Counter(mind[mind.minD<5].country).most_common()                                                                               
Out[72]: 
[('England', 198),
 ('USA', 142),
 ('China', 74),
 ('Scotland', 63),
 ('UAE', 57),
 ('Spain', 47),
 ('Iceland', 36),
 ('Australia', 31),
 ('Denmark', 29),
 ('Wales', 22),
 ('Russia', 18),
 ('Switzerland', 16),
 ('NewZealand', 15),
 ('Germany', 15),
 ('Belgium', 13),
 ('Sweden', 13),
 ('India', 12),
 ('Oman', 11),
 ('SouthAfrica', 10),
 ('Canada', 8),
 ('Indonesia', 4),
 ('Vietnam', 3),
 ('Italy', 3),
 ('Qatar', 2),
 ('Colombia', 2),
 ('Senegal', 2),
 ('Chile', 2),
 ('Austria', 1),
 ('Jordan', 1),
 ('SaudiArabia', 1),
 ('Turkey', 1),
 ('Israel', 1),
 ('Latvia', 1),
 ('Kazakhstan', 1),
 ('France', 1),
 ('MN908947', 1),
 ('Taiwan', 1)]

In [73]:                                                                                                                                
"""
