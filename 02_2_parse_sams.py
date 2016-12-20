import logging
from mirnylib import h5dict, genome
import argparse
import sys
from subprocess import call
from hiclib import mapping
from hiclib import fragmentHiC
import os
logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument("basename")
parser.add_argument("chunkNumber")
args = parser.parse_args()
basename = args.basename
chunk = args.chunkNumber
print(basename)

# B. Parse the mapped sequences into a Python data structure,
#    assign the ultra-sonic fragments to restriction fragments.
reads_file = '/exports/eddie/scratch/s1529682/processed/'+basename+'_'+chunk+'_mapped_reads.hdf5'
fragments_file = '/exports/eddie/scratch/s1529682/processed/'+basename+'_'+chunk+'_fragment_dataset.hdf5'
mapped_reads = h5dict.h5dict(reads_file)
genome_db    = genome.Genome('../genomes/mm9/fasta', readChrms=['#','X'])
def func():
    mapping.parse_sam(
        sam_basename1='/exports/eddie/scratch/s1529682/bams/'+basename+'_1.fq.gz'+chunk,
        sam_basename2='/exports/eddie/scratch/s1529682/bams/'+basename+'_2.fq.gz'+chunk,
        out_dict=mapped_reads,
        genome_db=genome_db, 
        enzyme_name='DpnII')
    fragments = fragmentHiC.HiCdataset(
        filename=fragments_file,
        genome=genome_db,
        maximumMoleculeLength=700,
        mode='w')
    
    # Load the parsed reads into the HiCdataset. The dangling-end filter is applied
    # at this stage, with maximumMoleculeLength specified at the initiation of the 
    # object.
    fragments.parseInputData(dictLike=reads_file)
         
try:
    func()
except ImportError:
    sys.stderr.write('Broken, restarting\n')
    call('ssh headnode1.ecdf.ed.ac.uk "cd /exports/igmm/datastore/wendy-lab/ilia/scripts; qsub 02_2_launch_parsing.sh ' + args.basename + '"', shell=True)
