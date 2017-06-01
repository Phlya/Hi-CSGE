# -*- coding: utf-8 -*-
"""
Created on Tue May 23 14:03:34 2017

@author: s1529682
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import cooler
import pandas as pd
from functools import partial
from multiprocessing import cpu_count, Pool
import itertools
import matplotlib.pyplot as plt
import seaborn as sns

def get_mids(intervals):
    intervals = intervals.sort_values(['Chrom', 'Start']).reset_index(drop=True)
    mids = np.round((intervals['End']+intervals['Start'])/2).astype(int)
    mids = pd.DataFrame({'Chrom':intervals['Chrom'], 'Mids':mids})
    return mids

def get_combinations(mids):
    combs = []
    for chrom in set(mids['Chrom']):
        current = np.array(list(itertools.combinations(mids[mids['Chrom']==chrom]['Mids'], 2)))
        if current.shape[0] < 2:
            continue
        combs.append(pd.DataFrame({'Chrom':chrom, 'Start':current[:,0], 'End':current[:,1]}))
    return pd.concat(combs).sort_values(['Chrom', 'Start', 'End']).dropna()

def controlLoops(df):
    """
    Creates "control" positions for loops

    :param df: loop (or domain) dataframe
    :return: ten copies of loop dataframe shifted randomly by 100kb to 1100kb
    """
    print('Making ctrls')
    dfs = []
    for i in range(10):
        df = df.copy()
        ran = (np.random.random(len(df)) + 0.1) * 1000000
        if i % 2 == 1:
            ran = -ran
        df["Start"] = df["Start"] + ran.astype(int)
        df["End"] = df["End"] + ran.astype(int)
        dfs.append(df)
    df = pd.concat(dfs)
    print('Made ctrl positions')
    return df

def get_vals(mychr, c, combinations):
    resolution = c.info['bin-size']
    vals = []
    current = combinations[combinations["Chrom"] == mychr]
    if not len(current) > 0:
        print('No data for %s' % mychr)
        return []
    print(mychr)
    data = c.matrix(sparse=True, balance=True).fetch(mychr).tocsr()
    for start, end in zip(current["Start"].values, current["End"].values):
        stBin, endBin = start//resolution, end//resolution
        try:
            vals.append([mychr, start, end, data[stBin, endBin]])
        except IndexError:
            continue
    return pd.DataFrame(vals, columns=['Chrom', 'Start', 'End', 'Val'])

def allVals(combinations, filename, pad = 8):
    c = cooler.Cooler(filename)
    f = partial(get_vals, c=c, combinations=combinations)
    p = Pool(cpu_count())
    myvals = p.map(f, c.chromnames[:-1])
    p.close()
    return pd.concat(myvals)
    
def allValsWithControl(combinations, filename):
    myvals =  allVals(combinations, filename)
    myctrls = allVals(controlLoops(combinations), filename)
    return myvals, myctrls
    
if __name__ == "__main__":
    import sys
    args = sys.argv
    baselist = args[1]
    coolfile = args[2]
    print(baselist)
    print(coolfile)
    bases = pd.read_csv(baselist, sep='\t',
                        names=['Chrom', 'Start', 'End'])
    combs = get_combinations(get_mids(bases))
    print('Now getting values...')
    loops = allVals(combs, coolfile)
    print('Got loops')
    ctrls = allVals(controlLoops(combs), coolfile)
    print('Got controls, saving')
    print(loops.shape)
    print(ctrls.shape)
    name = coolfile.split('/')[-1].split('_')[0]+' '+'_'.join(baselist.split('/')[-1].split('_')[:-1])
    loops.to_csv('../output/PcG_scaling/'+name+'.txt', sep='\t', index=False)
    ctrls.to_csv('../output/PcG_scaling/ctrl_'+name+'.txt', sep='\t', index=False)

    plt.loglog()
    plt.scatter((ctrls['End']-ctrls['Start'])//5000*5000, ctrls['Val'], label='Ctrls', c='r', alpha=0.1)
    plt.scatter((loops['End']-loops['Start'])//5000*5000, loops['Val'], label='PcG', c='b', alpha=0.1)
    plt.gca().set_aspect('equal')
    plt.legend()
    plt.title(name)
    plt.savefig('../output/PcG_scaling/%s_scaling.png' % name, dpi=300)
    ctrls['group'] = pd.cut((ctrls['End']-ctrls['Start'])//5000*5000, 10**np.arange(4, 9, 0.5))
    loops['group'] = pd.cut((loops['End']-loops['Start'])//5000*5000, 10**np.arange(4, 9, 0.5))
    ctrls['Kind'] = 'Ctrl'
    loops['Kind'] = 'CGI_PcG'
    d = pd.concat([loops, ctrls])
    d['Val'][d['Val']==0]=np.nan
    d = d.dropna()
    dists = d['group'].str.replace('(', '').str.replace(']', '').str.split(', ', expand=True).astype(float)
    d['dist'] = np.round(10**np.log10(dists).mean(axis=1), 1)
    plt.figure(figsize=(16, 16))
    plt.loglog()
    sns.pointplot(data=d, x='dist', y='Val', hue='Kind', estimator=np.mean)
    plt.gca().set_aspect('equal')
    plt.setp(plt.gca().get_xticklabels(), rotation=70)
    plt.savefig('../output/PcG_scaling/%s_scaling_pointplot.png' % name, dpi=300)
    plt.close()
