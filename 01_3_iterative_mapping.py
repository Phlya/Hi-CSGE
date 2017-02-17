#!/exports/applications/apps/SL7/anaconda/2.3.0/bin/python

import os
import sys
import logging
import argparse
from os import path
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("file", help="fastq file to map")
args = parser.parse_args()
sys.stderr.write(args.file+'\n')
logging.basicConfig(level=logging.DEBUG)

def func():
    #if not os.path.exists('tmp/'):
    #    os.mkdir('tmp/')
    
    # Map the reads iteratively.
    from hiclib import mapping
    #from mirnylib import h5dict, genome    
    
    mapping.iterative_mapping(
        bowtie_path='bowtie2',
        bowtie_index_path='../genomes/mm9/index/mm9',
        fastq_path=args.file,
        out_sam_path=path.join('/exports/eddie/scratch/s1529682/bams/',
                               path.split(args.file)[1] + '.bam'),
        min_seq_len=25,
        len_step=5,
        nthreads=4,
        #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
        temp_dir=tempfile.gettempdir(),  # optional, keep temporary files here
        bowtie_flags='--very-sensitive')

try:
    func()
except ImportError:
    sys.stderr.write('Broken, restarting\n')
    cwd = os.getcwd()
    from subprocess import call
    call('ssh headnode1.ecdf.ed.ac.uk "cd %s; qsub 01_2_launch_mapping.sh ' + args.file + '"' % cwd, shell=True)
