from Bio import SeqIO
import json
from collections import defaultdict, Counter
import pandas as pd
import pdb
import matplotlib.pyplot as plt

"""Determine whether transmissions happend domestically or internationally:
When parsing the phylogeny recursively, determine the origin of the ancestor. This information is passed on to the recursive calls.
Thus, when reaching a leave, it can be determined, whether at the time of infection of the leave, the lineage was already 
deemed domestic or not.

"""
datadir = '/Users/ahenschel/Github/ncov/auspice/'
#auspiceFile = 'covid19uae_Gisaid0525_b1_b2.json'
#auspiceFile = 'covid19uae_Gisaid0625_b1b2_k1k_countryMax300.json'
#auspiceFile = 'covid19uae_2020-06-25_b1b2N_28K.json'
#auspiceFile = 'covid19uae_2020-06-25_b1-3QCNR_countryMax100_k10.json'
auspiceFile = 'covid19uae_2020-06-25_b1-3QCNR_countryMax100_k25.json'
threshold = 50
weekWnd = 7/365.2422
day=1/365.2422

''' 
leaf = {'branch_attrs': {'labels': {'aa': 'N: A346T'},
  'mutations': {'N': ['A346T'], 'nuc': ['C18312T', 'G29309A']}},
 'name': 'Iran/KHGRC-1.1-IPI-8206/2020',
 'node_attrs': {'author': {'author': 'Sirous Zeinali et al',
   'value': 'Sirous Zeinali et al'},
  'country': {'confidence': {'Iran': 1.0},
   'entropy': -1.000088900581841e-12,
   'value': 'Iran'},
  'country_exposure': {'value': 'Iran'},
  'div': 0.00013499373540580567,
  'division': {'value': 'Semnan'},
  'division_exposure': {'value': 'Semnan'},
  'gisaid_epi_isl': {'value': 'EPI_ISL_442523'},
  'host': {'value': 'Human'},
  'num_date': {'confidence': [2020.1871584699454, 2020.1871584699454],
   'value': 2020.1871584699454},
  'originating_lab': {'value': 'Pasteur Institute of Iran'},
  'region': {'confidence': {'Asia': 0.9999999999999999},
   'entropy': -9.99866855976916e-13,
   'value': 'Asia'},
  'region_exposure': {'value': 'Asia'},
  'submitting_lab': {'value': 'Kawsar Human Genetic Research Company'},
  'url': 'https://www.gisaid.org'}}'''
"""auspice['tree']['node_attrs']

 {'country': {'confidence': {'Iran': 0.4357794187238073,
   'Jordan': 0.0526481309841954,
   'Kuwait': 0.297467544949543,
   'Turkey': 0.056334622591885826},
  'entropy': 1.7807770967169234,
  'value': 'Iran'},
 'div': 6.754749500105727e-05,
 'num_date': {'confidence': [2019.9277673878228, 2020.1147626262812],
  'value': 2020.0066895087193},
 'region': {'confidence': {'Asia': 0.9993530814775907},
  'entropy': 0.005446459183787267,
  'value': 'Asia'}}"""

mutationRec = []

def lineageMut(clade, inherMutations=[], protein='S', pos='614', origin='?', confidenceOrig=0):
    mutations = clade['branch_attrs'].get('mutations', {}).get(protein, [])
    location = clade['node_attrs']['country']['value']
    num_date = clade['node_attrs']['num_date']['value']
    confidence = clade['node_attrs']['country']['confidence'][location]
    transmissionType = 'global' if origin != location else 'local'
    accmutation = inherMutations + mutations #ut614
    mutPos = [mut for mut in accmutation if len(mut) == len(pos) + 2 and mut[1:-1] == pos]
    mut = mutPos[-1] if mutPos else None
    mutationRec.append((clade['name'], mut, num_date, location, confidence, transmissionType,
                        origin, confidenceOrig, 'children' in clade))
    #m614 = [mut for mut in mutations if pos in mut]
    #mut614 = [mut for mut in accmutation if len(mut) == len(pos) + 2 and mut[1:-1] == pos]
    #if mutations and m614:
    #    import pdb; pdb.set_trace()



    #if not 'children' in clade and (clade['name'].startswith('UAE') or clade['name'].startswith('UnitedArab')):
    #    print(clade['name'], accmutation)
    for child in clade.get('children', []):
        lineageMut(child, accmutation, protein, pos, location, confidence)


def lineage(clade, level=0, origin='?', confidenceOrig=0, countInternal=True):
    location = clade['node_attrs']['country']['value']
    num_date = clade['node_attrs']['num_date']['value']
    confidence = clade['node_attrs']['country']['confidence'][location]
    transmissionType = 'global' if origin != location else 'local'

    if countInternal or not 'children' in clade:
        transmissions.append((transmissionType, num_date, location, confidence, origin, confidenceOrig))

    for child in clade.get('children',[]):
        lineage(child, level+1, location, confidence)

def slidingWindow(table, windowSize=weekWnd*2, stride=weekWnd/2):
    table = table.sort_values(by='num_date')
    w1 = table['num_date'].iloc[0]
    w2 = w1 + windowSize
    tms = []
    while w2 < table['num_date'].iloc[-1]:
        glb = len(table[(w1 <= table.num_date) & (table.num_date <= w2) & (table.type=='global')])
        loc = len(table[(w1 <= table.num_date) & (table.num_date <= w2) & (table.type=='local')])
        mid = (w1 + w2) / 2
        tms.append((w1, w2, convertNumDate(mid), glb, loc, convertNumDate(w2)))
        w1 += stride
        w2 += stride
    return pd.DataFrame(tms, columns='ws1 ws2 mid global local date'.split())


def cumWindow(table, mutation, stride=day):
    table = table.sort_values(by='num_date')
    w1 = table['num_date'].iloc[0]
    w2 = w1 + stride
    tms = []
    while w2 < table['num_date'].iloc[-1]:
        tmp = table[table.num_date <= w2]
        total = len(tmp)
        mut = len(tmp[tmp.mutation==mutation])
        tms.append((w2, mut, total-mut, total,  mut/total, convertNumDate(w2)))
        w2 += stride
    return pd.DataFrame(tms, columns='num_date mut wt total ratio date'.split())

def convertNumDate(date): ##probably there is a function for this conversion somewhere...
    year = int(date)
    caldays = [31,28,31,30,31,30,31,31,30,31,30,31]
    if year%4==0 and year%100!=0: caldays[1] = 29
    days = (date - year) * sum(caldays)
    monthCount = 0
    while True:
        if days > caldays[monthCount]:
            days -= caldays[monthCount]
            monthCount += 1
        else: break
    return pd.Timestamp("%s-%02d-%02d" % (year, monthCount+1, int(days)+1))

def ratio(x, reportThr=-2):
    if x['global'] == x['local'] == 0 or max(x['global'], x['local'])< reportThr: return None
    if x['local'] == 0: return 1
    if x['global'] == 0: return 0
    return x['global']/(x['local'] + x['global'])

with open(f'{datadir}{auspiceFile}') as j:
    auspice = json.load(j)
transmissions = []
subclade = auspice['tree'] #['children'][0]['children'][1]
#lineage(subclade)

lineageMut(subclade)
mutation = 'D614G'
df = pd.DataFrame(mutationRec)
df.columns = ['name', 'mutation', 'num_date', 'location', 'confidence', 'trnsType', 'origin', 'origConf', 'internal']
#df1 = df[df.mutation==mutation]
for country, countryTable in df.groupby('location'):
    #if country != 'United Arab Emirates': continue
    if len(countryTable)<20: continue
    tms = cumWindow(countryTable, mutation)

    print(country, len(countryTable), len(tms))
    fig, ax = plt.subplots()
    ax.plot(tms['date'], tms['ratio'])
    ax.set_title(f'{country} {mutation} (acc.)')
    fig.show()
    fig.savefig(f'PlotsMut/mut_{mutation}_{country}.svg')
    fig.savefig(f'PlotsMut/mut_{mutation}_{country}.png')
    fig.clf()
    del fig



''' 
df = pd.DataFrame(transmissions, columns='type num_date location confidence origin conf_orig'.split())
df = df[(df.num_date<2020.5) & (df.num_date>2020.0)]


page = 0

for country, countryTable in df.groupby('location'):
    #if country != 'United Arab Emirates': continue

    if len(countryTable) < threshold: continue
    fig, ax = plt.subplots()
    print(country, len(countryTable))
    window = 14
    tms = slidingWindow(countryTable, window*day, stride=day)
    tms['ratio'] = tms.apply(ratio, axis=1)

    ax.plot(tms['date'], tms['global'], tms['date'], tms['local'])
    if country=='United Arab Emirates':
        ax.axvline('2020-03-25', color='k', linestyle='-.')
        ax.axvline('2020-03-18', color='k', linestyle=':')
    #ax2 = ax.twinx()
    #ax2.plot(tms['date'], tms['ratio'], 'g--')
    ax.set_title(country)
    plt.xticks(rotation=45)
    figname = f'Plots2/domVsInt2_{country}_{page:02d}'
    print(f'saving {figname}')
    fig.savefig(f'{figname}.svg',bbox_inches = "tight")
    fig.savefig(f'{figname}.png',bbox_inches = "tight")
    fig.clf()
#plt.show()
'''




