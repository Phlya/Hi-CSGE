# -*- coding: utf-8 -*-
import numpy as np
import cooler
import pandas as pd
import itertools
#from joblib import Parallel, delayed
from multiprocessing import Pool
from functools import partial
from scipy.misc import comb
import os


def get_mids(intervals):
    intervals = intervals.sort_values(['Chromosome', 'Start'])
    intervals = intervals[intervals['Chromosome'].isin(['chr'+ str(i) for i in list(range(1, 20))+['X']])].reset_index(drop=True)
    mids = np.round((intervals['End']+intervals['Start'])/2).astype(int)
    mids = pd.DataFrame({'Chromosome':intervals['Chromosome'], 'Mids':mids}).drop_duplicates()
    return mids

def get_combinations(mids, res, local=False):
    m = (mids['Mids']//res).astype(int).values
    if local:
        for i in m:
            yield i, i
    else:
        for i in itertools.combinations(m, 2):
#            if abs(i[1]-i[0])<10**6/res:#############################################
                yield i

def controlLoops(midcombs, res, minshift=10**5, maxshift=10**6, nshifts=1):
    minbin = minshift//res
    maxbin = maxshift//res
    for start, end in midcombs:
        for i in range(nshifts):
            shift = np.random.randint(minbin, maxbin)*np.sign(np.random.random()-0.5).astype(int)
            yield start+shift, end+shift

def averageLoops(chrom, mids, c, pad = 7, ctrl=False, local=False):
    mymap = np.zeros((2*pad + 1, 2*pad + 1), np.float64)
    data = c.matrix(sparse=True, balance=True).fetch(chrom).tocsr()
    
    current = mids[mids["Chromosome"] == chrom]
    if not len(current) > 1:
        mymap.fill(np.nan)
        return mymap
    if ctrl:
        current = controlLoops(get_combinations(current, c.binsize, local), c.binsize)
    else:
        current = get_combinations(current, c.binsize, local)
    n = 0
    for stBin, endBin in current:
        if not local and abs(endBin - stBin) < pad*2:
            continue
        try:
            mymap += np.nan_to_num(data[stBin - pad:stBin + pad+1, endBin - pad:endBin + pad+1].toarray())
            n += 1
        except (IndexError, ValueError) as e:
            continue
    return mymap/n

def averageLoopsWithControl(mids, filename, pad, nproc, chroms, local):
    p = Pool(nproc)
    c = cooler.Cooler(filename)
    #Loops
    f = partial(averageLoops, mids=mids, c=c, pad=pad, local=local)
    loops = map(f, chroms)
    #Controls
    f = partial(averageLoops, mids=mids, c=c, pad=pad, ctrl=True, local=local)
    ctrls = map(f, chroms)
    p.close()
    return loops, ctrls
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("coolfile", type=str)
    parser.add_argument("baselist", type=str)
    parser.add_argument("--pad", default=7, type=int, required=False)
    parser.add_argument("--n_proc", default=0, type=int, required=False)
    parser.add_argument("--excl_chrs", default='chrY,chrM', type=str, required=False)
    parser.add_argument("--incl_chrs", default='all', type=str, required=False)
    parser.add_argument("--local", action='store_true', default=False, required=False)
    parser.add_argument("--outdir", default='../output/new pileups 5Kb', type=str, required=False)
    args = parser.parse_args()
    print(args.coolfile)
    print(args.baselist)
    if args.n_proc==0:
        nproc=-1
    else:
        nproc=args.n_proc
    
    c = cooler.Cooler(args.coolfile)
    if args.incl_chrs=='all':
        incl_chrs = c.chromnames
    else:
        incl_chrs = args.incl_chrs.split(',')
    chroms = cooler.Cooler(args.coolfile).chromnames
    fchroms = []
    for chrom in chroms:
        if chrom not in args.excl_chrs.split(',') and chrom in incl_chrs:
            fchroms.append(chrom)
            
    bases = pd.read_csv(args.baselist, sep='\t',
                        names=['Chromosome', 'Start', 'End'])
    mids = get_mids(bases)
    loops, ctrls = averageLoopsWithControl(mids, args.coolfile, args.pad, nproc, fchroms, args.local)
    w = mids.groupby('Chromosome').size()
    w = w.loc[fchroms].values
    w = comb(w, 2)
    lp = np.average([lp for lp in loops if not np.all(np.isnan(lp))], weights=w, axis=0)
    ctrl = np.average([ctrl for ctrl in ctrls if not np.all(np.isnan(ctrl))], weights=w, axis=0)
    loop = lp/ctrl
    np.savetxt(os.path.join(args.outdir, args.coolfile.split('/')[-1].split('_')[0]+'_'+'_'.join(args.baselist.split('/')[-1].split('_')[:-1])+'.np.txt'), loop)
