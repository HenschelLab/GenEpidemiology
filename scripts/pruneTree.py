from ete3 import Tree
import pandas as pd
from Bio import SeqIO
import sys
keepCountries = ['UAE']
# run ../scripts/pruneTree.py ../results/tree_raw_2020-06-25_b1b2N.nwk 

#t = Tree("tree_closest300.nwk", format=1) 
t = Tree(sys.argv[-1], format=1) 
chineseCities = [ 'Anhui',  'Beijing', 'Chongqing', 'Foshan', 'Fujian', 'Fuyang', 'Ganzhou','Guangdong', 'Guangzhou', 'Hangzhou', 'Hefei', 'Henan', 'HongKong', 
                  'Jian', 'Jiangsu', 'Jiangxi', 'Jingzhou', 'Jiujiang',  'Lishui',  'MN908947_S',  'NanChang', 'Nanchang',
                  'Pingxiang',  'Shandong', 'Shaoxing', 'Shanghai', 'Shangrao', 'Shenzhen', 'Sichuan',   'Tianmen', 'Wuhan', 'Wuhan-Hu-1', 'Xinyu', 'Yunnan', 'Zhejiang']

                  #'canine', 'cat', 'env', 'hCoV-19_S', 'mink', 'tiger'

#colors = [cc.cm.CET_D1A(c) for c in np.linspace(0,1.05, maxScore)] ##WTF goes wrong here: NoneType error
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
        return year in [2019,2020] and 1<=month<=12 and 1 <= day <= 31
    except:
        return False


fastafile = "../data/sequences_2020-06-25_b1b2.fasta"

df = pd.DataFrame([seqStats(rec) for rec in SeqIO.parse(fastafile, "fasta")])
df.columns = 'descr country length gaps N Nclipped'.split()

nonhuman = ['canine', 'cat', 'env', 'hCoV-19_S', 'mink', 'tiger', 'pangolin', 'bat', 'mouse'] 
df1 = df[(30000>df.length) & (df.length>29000) & (df.Nclipped < 30) & (~df.country.isin(nonhuman))]  

leaves = t.get_leaf_names() 
retain = set(leaves).intersection(df1.descr)


md = pd.read_csv("../data/metadata_2020-06-25b1b2.tsv", sep='\t', parse_dates=['date'])
md['datecheck'] = md.apply(lambda x:checkDate(x['date']), axis=1)
mdQ = md[md.datecheck]
mdQ = mdQ.set_index('strain')
highQual = set(mdQ.index.values).intersection(leaves)
mdQ = mdQ.loc[highQual]
#mdQ.to_csv('../data/metadata_2020-06-25b1b2N.tsv', sep='\t')

#retain2 = retain.intersection(mdQ.index.values)
#t.prune(retain2)
#t.write(format=1, outfile="../results/tree_pruned_2020-06-25_b1b2Na.nwk") 

