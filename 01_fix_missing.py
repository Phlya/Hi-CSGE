#!/exports/applications/apps/SL7/anaconda/2.3.0/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 10:33:33 2016

@author: s1529682
"""

import glob
import subprocess

bams = '/exports/eddie/scratch/s1529682/bams/*'
fastq = '/exports/eddie/scratch/s1529682/fastq/split/*'
ffiles = glob.glob(fastq)
bfiles = glob.glob(bams)
j = 0
for ff in ffiles:
    bn = ff.split('/')[-1]
    sample = bn.split('.')[0]
    if sample.split('_')[0][-1]=='t':
        n = 11
    else:
        n = 16
    bfs = [i for i in bfiles if i.split('/')[-1].startswith(bn)]
    if len(bfs)<n:
#        subprocess.call('ssh headnode1.ecdf.ed.ac.uk "cd /exports/igmm/datastore/wendy-lab/ilia/scripts; qsub 01_2_launch_mapping.sh %s"' % ff, shell=True)
        print ff, len(bfs)
