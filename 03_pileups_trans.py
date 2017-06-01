# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import cooler
import pandas as pd
from functools import partial
import itertools
from multiprocessing import cpu_count, Pool
from joblib import delayed, Parallel

def get_mids(intervals):
    intervals = intervals.sort_values(['Chrom', 'Start'])
    intervals = intervals[intervals['Chrom'].isin(['chr'+ str(i) for i in list(range(1, 20))+['X']])].reset_index(drop=True)
    mids = np.round((intervals['End']+intervals['Start'])/2).astype(int)
    mids = pd.DataFrame({'Chrom':intervals['Chrom'], 'Mids':mids})
    return mids

def get_trans_combinations(mids):
    combs = []
    for chrom1, chrom2 in itertools.combinations(set(mids['Chrom']), 2):
        current1 = np.array(mids[mids['Chrom']==chrom1]['Mids'])
        current2 = np.array(mids[mids['Chrom']==chrom2]['Mids'])
        if current1.shape[0] < 1 or current2.shape[0] < 1:
            continue
        prod = np.array(list(itertools.product(current1, current2)))
        combs.append(pd.DataFrame({'Chrom1':chrom1, 'Chrom2':chrom2, 'Start':prod[:,0], 'End':prod[:,1]}))
    return pd.concat(combs).sort_values(['Chrom1', 'Chrom2', 'Start', 'End']).dropna()

def controlLoops(df):
    """
    Creates "control" positions for loops

    :param df: loop (or domain) dataframe
    :return: ten copies of loop dataframe shifted randomly by 100kb to 1100kb
    """
    dfs = []
    for i in range(1):
        df = df.copy()
        ran = (np.random.random(len(df)) + 0.1) * 1000000
        if i % 2 == 1:
            ran = -ran
        df["Start"] = df["Start"] + ran
        ran = (np.random.random(len(df)) + 0.1) * 1000000
        if i % 2 == 1:
            ran = -ran
        df["End"] = df["End"] + ran
        dfs.append(df)
    df = pd.concat(dfs)
    return df

def get_submaps(mychrs, c, mids, pad, ctrl=False):
    mychr1, mychr2 = mychrs
    print(mychr1, mychr2)
    mymap = np.zeros((2 * pad, 2 * pad), np.float64)
    resolution = c.info['bin-size']
    data = c.matrix(sparse=True, balance=True).fetch(mychr1, mychr2).tocsr()
    data[data!=data] = 0
    data /= data.mean()
    current = mids[(mids["Chrom"] == mychr1) |
                   (mids["Chrom"] == mychr2)]
    current = get_trans_combinations(current).dropna()
    if ctrl:
        current = controlLoops(current)
    if not len(current) > 0:
        return mymap
    n = 0
    for stBin, endBin in zip(current["Start"].values, current["End"].values):
        stBin  = int(stBin//resolution)
        endBin  = int(endBin//resolution)
        if stBin - pad < 0 or stBin + pad > data.shape[0]:
            continue
        if endBin + pad > data.shape[1] or endBin - pad < 0:
            continue
        mymap = mymap + np.nan_to_num(data[stBin - pad:stBin + pad, endBin - pad:endBin + pad].toarray())
        n += 1
    return mymap/n

def averageLoops(mids, filename, pad = 8):
    c = cooler.Cooler(filename)
    f = partial(get_submaps, c=c, mids=mids, pad=pad)
#    p = Pool(cpu_count())
#    mymaps = p.map(f, itertools.combinations(c.chromnames[:-1], 2))
    mymaps = Parallel(-1)(delayed(f)(i) for i in itertools.combinations(c.chromnames[:-1], 2))
#    f = partial(get_submaps, c=c, mids=mids, pad=pad, ctrl=True)
#    myctrls = Parallel(-1)(delayed(f)(i) for i in itertools.combinations(c.chromnames[:-1], 2))
    return mymaps#, myctrls

if __name__ == "__main__":
    import sys
    args = sys.argv
    baselist = args[1]
    coolfile = args[2]
    print(baselist)
    print(coolfile)
    try:
        pad = int(args[3])
    except IndexError:
        pad = 8

    bases = pd.read_csv(baselist, sep='\t',
                        names=['Chrom', 'Start', 'End'])
    print('Read bed file')
    mids = get_mids(bases)
    print('Made mids')
    loops = averageLoops(mids, coolfile, pad)
    print('Made loops, saving')
#    loop = np.array([lp/ctrl for lp, ctrl in zip(loops, ctrls)]).mean(axis=0)
    loop = np.array(loops).mean(axis=0)
    r = int(cooler.Cooler(coolfile).info['bin-size']//1000)
    np.savetxt('../output/pileups %sKb/trans/'%r+coolfile.split('/')[-1].split('_')[0]+'_'+'_'.join(baselist.split('/')[-1].split('_')[:-1])+'.np.txt', loop)
    print('All done!')
