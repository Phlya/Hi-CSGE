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
    bfs = [i for i in bfiles if i.split('/')[-1].startswith(bn)]
    if len(bfs)<16:
        subprocess.call('qsub 01_2_launch_mapping.sh %s' % ff, shell=True)
#        j += 1
#        print ff 
