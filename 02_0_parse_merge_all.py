# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 18:20:13 2016

@author: s1529682
"""
import sys
import os
import subprocess

path = '/exports/eddie/scratch/s1529682/bams'
add_basenames = ['serum', '2i', 'R2i', 'KO', 'I53A']
files = os.listdir(path)
basenames = set([f.split('_')[0] for f in files])
#for bn in basenames:
#    subprocess.call('qsub launch_parsing.sh ' + bn, shell=True)
chunks = {}
for bn in basenames:
    names = [f for f in files if f.startswith(bn)]
    chunks[bn] = set([name.split('gz')[1].split('.')[0] for name in names])
    
for bn in chunks:
    for chunk in chunks[bn]:
        subprocess.call('qsub -N parsing%s 02_1_launch_parsing.sh %s %s' % (bn, bn, chunk), shell=True)
    subprocess.call('qsub -hold_jid parsing%s -N merging%s 02_3_launch_merging.sh %s'  % (bn, bn, bn), shell=True)

dependencies = {}
for addbn in add_basenames:
    dependencies[addbn] = []
    for bn in chunks:
        if bn.startswith(addbn):
            dependencies[addbn].append(bn)

for addbn in add_basenames:
    subprocess.call('qsub -hold_jid %s -N merging%s 02_3_launch_merging.sh %s' % (','.join(dependencies[addbn]), addbn, addbn), shell=True)
