## removing few leaves from timetree using ETE
## use e.g. conda environment nextstrain
import os
from ete3 import Tree
t = "../results/tree_2020-06-25_b1-3QCNR_countryMax100_k25.nwk" 
out = "/tmp/prunedTree.nwk"

t1 = Tree(t, format=1)
prun=['Turkey/6224-Ankara1034/2020', 'Canada/ON_MU-S53/2020', 'SouthAfrica/R02606/2020',
      'BosniaandHerzegovina/03_Tuzla/2020'] 
leaves = [l for l in t1.iter_leaf_names()]
l1 = set(leaves).difference(prun)
t1.prune(l1) 
t1.write(format=1, outfile=out)

os.system(f"mv {t} {t}_bak")
os.system(f"mv {out} {t}")
