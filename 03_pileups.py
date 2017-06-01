# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import cooler
#from mirnylib.genome import Genome
#from hiclib import hicShared
import pandas as pd
from mirnylib.numutils import coarsegrain
import itertools
from joblib import Parallel, delayed
from functools import partial
#mm9 = Genome('~/Documents/PhD/Hi-C/genomes/mm9/fasta/', readChrms=["#", "X"])

def get_mids(intervals):
    intervals = intervals.sort_values(['Chromosome', 'Start'])
    intervals = intervals[intervals['Chromosome'].isin(['chr'+ str(i) for i in list(range(1, 20))+['X']])].reset_index(drop=True)
    mids = np.round((intervals['End']+intervals['Start'])/2).astype(int)
    mids = pd.DataFrame({'Chromosome':intervals['Chromosome'], 'Mids':mids})
    return mids

def get_combinations(mids):
    combs = []
    for chrom in set(mids['Chromosome']):
        current = np.array(list(itertools.combinations(mids[mids['Chromosome']==chrom]['Mids'], 2)))
        combs.append(pd.DataFrame({'Chromosome':chrom, 'Start':current[:,0], 'End':current[:,1]}))
    return pd.concat(combs).sort_values(['Chromosome', 'Start', 'End'])

def controlLoops(df):
    """
    Creates "control" positions for loops

    :param df: loop (or domain) dataframe
    :return: ten copies of loop dataframe shifted randomly by 100kb to 1100kb
    """
    dfs = []
    for i in range(10):
        df = df.copy()
        ran = (np.random.random(len(df)) + 0.1) * 1000000
        if i % 2 == 1:
            ran = -ran
        df["Start"] = df["Start"] + ran
        df["End"] = df["End"] + ran
        dfs.append(df)
    df = pd.concat(dfs)
    return df

def averageLoops(loopPositions, filename, pad = 8):
    c = cooler.Cooler(filename)
#    mygen = mm9
    resolution = c.info['bin-size']

    mymaps = []

    for mychr in c.chromnames[:-1]:
        mymap = np.zeros((2 * pad, 2 * pad), np.float64)
        print(mychr)

        #data = myd.get_dataset("{0} {0}".format(mychr))
        data = c.matrix(sparse=True, balance=True).fetch(mychr).tocsr()
        # plt.imshow(data[:2000, :2000])
        # plt.show()
        current = loopPositions[loopPositions["Chromosome"] == mychr]
        if not len(current) > 0:
            continue
        current = np.floor(current[['Start', 'End']]/resolution).astype(int)

        for stBin, endBin in zip(current["Start"].values, current["End"].values):
            if abs(stBin - endBin) < pad + 2:
                continue
            # print (stBin, endBin)
            if stBin - pad < 0:
                continue
            if endBin + pad > data.shape[0]:
                continue
            # plt.imshow(data[stBin - pad:stBin + pad, endBin - pad:endBin + pad])
            # plt.show()
            mymap = mymap + np.nan_to_num(data[stBin - pad:stBin + pad, endBin - pad:endBin + pad].toarray())
        mymaps.append(mymap)
    for i in mymaps:
        assert i.shape == (2 * pad, 2*pad)
    return mymaps

def get_submap(start, end, data, pad):
    return np.nan_to_num(data[start - pad:start + pad, end - pad:end + pad].toarray())

def averageLoopsParallel(loopPositions, filename, pad = 8):
    c = cooler.Cooler(filename)
#    mygen = mm9
    resolution = c.info['bin-size']

    mymaps = []

    for mychr in c.chromnames[:-1]:
#        mymap = np.zeros((2 * pad, 2 * pad), np.float64)
        print(mychr)

        #data = myd.get_dataset("{0} {0}".format(mychr))
        data = c.matrix(sparse=True, balance=True).fetch(mychr).tocsr()
        # plt.imshow(data[:2000, :2000])
        # plt.show()
        current = loopPositions[loopPositions["Chromosome"] == mychr]
        assert len(current) > 0
        current = np.floor(current[['Start', 'End']]/resolution).astype(int)
        
        bins = []
        for stBin, endBin in zip(current["Start"].values, current["End"].values):
            if abs(stBin - endBin) < 10:
                continue
            # print (stBin, endBin)
            if stBin - pad - 2 < 0:
                continue
            if endBin + pad -2 > data.shape[0]:
                continue
            bins.append((stBin, endBin))
            # plt.imshow(data[stBin - pad:stBin + pad, endBin - pad:endBin + pad])
            # plt.show()
        get = partial(get_submap, data=data, pad=pad)
        mymap = np.sum(Parallel(n_jobs=3)(delayed(get)(start, end) for start, end in bins), axis=0)
#        mymap = mymap + np.nan_to_num(data[stBin - pad:stBin + pad, endBin - pad:endBin + pad].toarray())
        mymaps.append(mymap)
    for i in mymaps:
        assert i.shape == (2 * pad, 2*pad)
    return mymaps

def averageLoopsWithControl(loopPositions, filename, cg=1, pad=8):
    mymaps =  averageLoops(loopPositions, filename, pad = pad)
    mymaps2 = averageLoops(controlLoops(loopPositions), filename, pad = pad)
    if cg != 1:
        mymaps = [coarsegrain(1. * i, cg) / (cg ** 2) for i in mymaps]
        mymaps2 = [coarsegrain(1. * i, cg) / (cg ** 2) for i in mymaps2]
    return mymaps, mymaps2
    
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
    
    try:
        cg = int(args[4])
    except IndexError:
        cg = 1
        
   
    bases = pd.read_csv(baselist, sep='\t',
                        names=['Chromosome', 'Start', 'End'])
    combs = get_combinations(get_mids(bases))
    loops, ctrls = averageLoopsWithControl(combs, coolfile, cg, pad)
    loop = np.array([lp/ctrl*10 for lp, ctrl in zip(loops, ctrls)]).mean(axis=0)
    np.savetxt('../output/'+coolfile.split('/')[-1].split('_')[0]+'_'+'_'.join(baselist.split('/')[-1].split('_')[:-1])+'.np.txt', loop)
