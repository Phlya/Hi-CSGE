import numpy as np
import pandas as pd
from peakdetect import peakdetect #https://gist.github.com/sixtenbe/1178136#file-peakdetect-py

def insul_diamond(A, extent=200):
    N = A.shape[0]
    score = np.zeros(N)
    for i in range(0, N):
        lo = max(0, i-extent)
        hi = min(i+extent, N)
        score[i] = np.nansum(A[lo:i, i:hi].toarray())
    score /= score.mean()
    return score

def up_down_contacts(A, extent=200):
    N = A.shape[0]
#     A = np.triu(A, 0)
    score = np.zeros(N)
    for i in range(0, N):
        lo = max(0, i-extent)
        hi = min(i+extent, N)
        sub = np.triu(A[lo:hi, lo:hi].toarray(), 2)
        score[i] = np.nansum(sub[0:extent, 0:extent])+np.nansum(sub[extent:-1, extent:-1])
    score /= score.mean()
    return score

def insul_score(A, extent=200):
    score = up_down_contacts(A, extent)/insul_diamond(A, extent)
    score[np.isinf(score)]=np.nan
    return score

def call_TADs(cool, name, func=insul_score, valleys=False, il=20, la=5, d=0.25, color='0,0,0'):
    results = []
    for chrom in list(range(1, 20))+['X']:
        print(chrom)
        data = cool.matrix(sparse=True).fetch('chr'+str(chrom)).tocsr()
        S = insul_score(data, il)
        S /= np.nanmean(S)
        peaks = peakdetect(S, lookahead=la, delta=d)[int(valleys)]
        peaks = np.array(peaks).T[0].astype(int)
        tads = pd.DataFrame(np.array(list(zip(peaks[:-1], peaks[1:]))), columns=['x1', 'x2'])*serum.binsize
        tads['chr1'] = chrom
        tads_annot = pd.concat([tads, tads.rename(columns={'chr1':'chr2', 'x1':'y1', 'x2':'y2'})], axis=1)['chr1 x1 x2 chr2 y1 y2'.split()]
        results.append(tads_annot)
    TADs_annot = pd.concat(results)
    TADs_annot['color'] = color
    TADs_annot['comment'] = '%s insul_score %s %s %s' % (name, il, la, d)
    return TADs_annot

def call_TADs_chrom(chrom, cool, func, valleys, il, la, d):
    print(chrom)
    data = cool.matrix(sparse=True).fetch(chrom).tocsr()
    S = func(data, il)
    S /= np.nanmean(S)
    peaks = peakdetect(S, lookahead=la, delta=d)[int(valleys)]
    peaks = np.array(peaks).T[0].astype(int)
    tads = pd.DataFrame(np.array(list(zip(peaks[:-1], peaks[1:]))), columns=['x1', 'x2'])*serum.binsize
    tads['chr1'] = chrom
    tads_annot = pd.concat([tads, tads.rename(columns={'chr1':'chr2', 'x1':'y1', 'x2':'y2'})], axis=1)['chr1 x1 x2 chr2 y1 y2'.split()]
    return tads_annot

def pcall_TADs(cool, name, func=insul_score, n_cores=-1, valleys=False, il=20, la=5, d=0.25, color='0,0,0'):
    #Good params for mouse genome at 10Kb resolution, takes around 1 min on intermediate computer with 4 cores
    f = partial(call_TADs_chrom, cool=cool, func=func, valleys=valleys, il=il, la=la, d=d)
    results = Parallel(n_cores)(delayed(f)(chrom) for chrom in cool.chromnames[:-1])
    TADs_annot = pd.concat(results)
    TADs_annot['color'] = color
    TADs_annot['comment'] = '%s insul_score %s %s %s' % (name, il, la, d)
    return TADs_annot
