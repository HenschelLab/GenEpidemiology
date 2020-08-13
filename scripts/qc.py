import pandas as pd
from Bio import SeqIO

chineseCities = [ 'Anhui',  'Beijing', 'Chongqing', 'Foshan', 'Fujian', 'Fuyang', 'Ganzhou','Guangdong', 'Guangzhou', 'Hangzhou', 'Hefei', 'Henan', 'HongKong', 
                  'Jian', 'Jiangsu', 'Jiangxi', 'Jingzhou', 'Jiujiang',  'Lishui',  'MN908947_S',  'NanChang', 'Nanchang',
                  'Pingxiang',  'Shandong', 'Shaoxing', 'Shanghai', 'Shangrao', 'Shenzhen', 'Sichuan',   'Tianmen', 'Wuhan', 'Wuhan-Hu-1', 'Xinyu', 'Yunnan', 'Zhejiang']

                  #'canine', 'cat', 'env', 'hCoV-19_S', 'mink', 'tiger'

def selectSequences(fromFile="../results/aligned.fasta", toFile="../results/sel_Nclipped0_10K.fasta", selCrit=lambda i,rec:True):
    counter=0
    with open(toFile, 'w') as w:
        for i, rec in enumerate(SeqIO.parse(fromFile, "fasta")):
            if selCrit(i, rec): #rec.description in selAll:#selLocal1000: # or rec.description.split('/')[0] in selectedLocations:
                if rec.id == 'Wuhan/Hu-1/2019': rec.id = 'MN908947'
                print(f">{rec.description}\n{rec.seq}", file=w)
                counter += 1
    print(f"Wrote {counter} sequences to {toFile}")

def getcountry(loc):
    if loc in chineseCities: return 'China'
    elif loc == 'UnitedArabEmirates': return 'UAE'
    elif loc == 'ITALY': return 'Italy'
    else: return loc
def seqStats(rec):
    return (rec.description, getcountry(rec.description.split('/')[0]),
            len(rec.seq), rec.seq.count('-'), rec.seq.count('N'), rec.seq[65:-73].count('N'))
def checkDate(date):
    try:
        year, month, day = list(map(int, date.split('-')))
        return year in [2019,2020] and 1<=month<=12 and 1 <= day <= 31 and (2019,12) <= (year, month) < (2020, 8) 
    except:
        return False

spec = "2020-06-25_b1-3"
#spec = "2020-06-25_b1b2"

fastafile = f"../results/aligned_{spec}.fasta"
fastafileQC = f"../results/aligned_{spec}QC.fasta"


df = pd.DataFrame([seqStats(rec) for rec in SeqIO.parse(fastafile, "fasta")])
df.columns = 'strain country length gaps N Nclipped'.split()

nonhuman = ['canine', 'cat', 'env', 'hCoV-19_S', 'mink', 'tiger', 'pangolin', 'bat', 'mouse'] 
df1 = df[(30000>df.length) & (df.length>29000) & (df.Nclipped < 30) & (~df.country.isin(nonhuman))]  

#leaves = t.get_leaf_names() 
#retain = set(leaves).intersection(df1.descr)


md = pd.read_csv(f"../data/metadata_{spec}.tsv", sep='\t') #, parse_dates=['date'])
md['datecheck'] = md.apply(lambda x:checkDate(x['date']), axis=1)
mdQ = md[md.datecheck]
#mdQ = mdQ.set_index('strain')
#highQual = set(mdQ.index.values).intersection(leaves)
#mdQ = mdQ.loc[highQual]
qc = pd.merge(mdQ, df1, on='strain', how='inner') 
qc[qc.columns[:-6]].to_csv(f'metadata_{spec}QC.tsv', sep='\t', index=False)  
selection = set(qc['strain'])
selCrit = lambda i, rec: rec.id in selection
selectSequences(fastafile, fastafileQC, selCrit=selCrit)
### NOTE!!!!
# still manually fixing Wuhan/Hu-1 in metadata and fasta file, TO BE Fixed!



