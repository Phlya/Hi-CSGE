from mirnylib import genome
from hiclib import fragmentHiC 
from glob import glob
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("sample")
args = parser.parse_args()
sample = args.sample

genome_db = genome.Genome('../genomes/mm9/fasta', readChrms=['#','X'])
genome_db.setEnzyme('DpnII')
path = '/exports/eddie/scratch/s1529682/processed/'

files = glob(path+sample+'*fragment_dataset.hdf5')

print(sample)
fragments = fragmentHiC.HiCdataset(
    filename=path+'merged/'+sample+'_fragment_dataset.hdf5',
    genome=genome_db,
    maximumMoleculeLength=750,
    mode='w',
    tmpFolder='/exports/eddie/scratch/s1529682/tmp')
fragments.merge(files)
print('Merged', sample)
fragments = fragmentHiC.HiCdataset(
    filename=path+'merged/'+sample+'_fragment_dataset_filtered.hdf5',
    genome=genome_db,
    maximumMoleculeLength=750,
    mode='w',
    tmpFolder='/exports/eddie/scratch/s1529682/tmp')
#fragments._sortData()
fragments.load(path+'merged/'+sample+'_fragment_dataset.hdf5')
fragments.filterDuplicates(mode='ram')
fragments.filterLarge(10000,10)   # DpnII
fragments.filterExtreme(cutH=0.0001, cutL=0)
print('Filtered', sample)

#fragments.filterRsiteStart(offset=5)
#fragments.filterLarge(cutsmall=50)
#fragments.writeFilteringStats()
#fragments.filterDuplicates(mode='hdd', tmpDir='../tmp')
#fragments.filterExtreme(cutH=0.005, cutL=0)

fragments.writeFilteringStats()
fragments.printMetadata(sample+'_stats.txt')
x, y = fragments.plotScaling(plot=False)
np.savetxt(sample + '_scaling.np.txt', np.vstack([x, y]))

#def make_cooler(frags, bins, cool_path):
#    bins = pd.read_csv(
#        bins,
#        sep='\t',
#        names=['chrom', 'start', 'end'],
#        dtype={'chrom': str})

    # Chrom sizes from bin table
#    chromtable = (
#        bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
#            .reset_index(drop=True)
#            .rename(columns={'chrom': 'name', 'end': 'length'})
#    )
#    chroms, lengths = list(chromtable['name']), list(chromtable['length'])
#    chromsizes = pd.Series(index=chroms, data=lengths)

    # Aggregate the contacts
#    chunksize = int(10e6)
#    with h5py.File(frags, 'r') as h5pairs, \
#         h5py.File(cool_path, 'w') as h5:
#        reader = cooler.io.HDF5Aggregator(h5pairs, chromsizes, bins, chunksize)
#        cooler.io.create(h5, chroms, lengths, bins, reader) # metadata, assembly)
#        pool = Pool(16)
#        bias = cooler.ice.iterative_correction(h5, map=pool.map)
#        h5opts = dict(compression='gzip', compression_opts=6)
#        h5['bins'].create_dataset('weight', data=bias, **h5opts)
#        h5.flush()
#        h5.close()

#for i in (5, 10, 25, 100, 1000)[::-1]:
#    i = str(i)
#    make_cooler(path+'merged/'+sample+'_fragment_dataset.hdf5', i+'Kb_bins_mm9.bed', path+'merged/'+sample+'_'+i+'Kb.cool')

#juicepath = path + 'merged/' + 'for_juicebox/' + sample + '_for_juicebox.txt'

#with open(juicepath, 'w') as f:
#    for str1, str2, chr1, chr2, pos1, pos2, dist1, dist2, mapq in zip(fragments.h5dict['strands1'], fragments.h5dict['strands2'],
#                                                                            fragments.h5dict['chrms1'], fragments.h5dict['chrms2'],
#                                                                            fragments.h5dict['cuts1'], fragments.h5dict['cuts2'],
#                                                                           fragments._getVector('dists1'), fragments._getVector('dists2'),
#                                                                           repeat(32)):
#	    chr1 = 'chr'+str(chr1+1) if chr1<19 else 'chrX'
#	    chr2 = 'chr'+str(chr2+1) if chr2<19 else 'chrX'
#	    frag1 = pos1 + dist1 * (str1 * 2 - 1) 
#	    frag2 = pos2 + dist1 * (str2 * 2 - 1) 
#	    str1 = int(str1) - 1
#	    str2 = int(str2) - 1
#	    f.write('\t'.join([str(i) for i in (str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, int(mapq), int(mapq))] + '\n'))
