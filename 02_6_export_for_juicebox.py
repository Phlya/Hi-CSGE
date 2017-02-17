from mirnylib import genome
from hiclib import fragmentHiC
from itertools import repeat
import argparse
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("sample")
args = parser.parse_args()
sample = args.sample

genome_db = genome.Genome('../genomes/mm9/fasta', readChrms=['#','X'])
genome_db.setEnzyme('DpnII')

path = '/exports/eddie/scratch/s1529682/processed/'
fragments = fragmentHiC.HiCdataset(
    filename='bla',
    inMemory=True,
    genome=genome_db,
    maximumMoleculeLength=750,
    mode='w',
    tmpFolder=tempfile.gettempdir())
fragments.load(path+'merged/'+sample+'_fragment_dataset_filtered.hdf5')

juicepath = path + 'merged/for_juicebox/' + sample + '_for_juicebox.txt'

with open(juicepath, 'w') as f:
    for str1, str2, chr1, chr2, pos1, pos2, dist1, dist2, mapq in zip(fragments.h5dict['strands1'], fragments.h5dict['strands2'],
                                                                            fragments.h5dict['chrms1'], fragments.h5dict['chrms2'],
                                                                            fragments.h5dict['cuts1'], fragments.h5dict['cuts2'],
                                                                           fragments._getVector('dists1'), fragments._getVector('dists2'),
                                                                           repeat(32)):
	    chr1 = 'chr'+str(chr1+1) if chr1<19 else 'chrX'
	    chr2 = 'chr'+str(chr2+1) if chr2<19 else 'chrX'
	    frag1 = pos1 + dist1 * (str1 * 2 - 1) 
	    frag2 = pos2 + dist1 * (str2 * 2 - 1) 
	    str1 = int(str1) - 1
	    str2 = int(str2) - 1
	    f.write('\t'.join([str(i) for i in (str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, int(mapq), int(mapq))]) + '\n')
