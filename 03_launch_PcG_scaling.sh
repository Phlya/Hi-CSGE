#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd
#$ -N PcG_scaling
#$ -l h_rt=06:00:00
#$ -l h_vmem=100G
#$ -pe sharedmem 4
. /etc/profile.d/modules.sh

module load anaconda/2.3.0
module load igmm/apps/hdf5/1.8.16
source activate Hi-C

python3 03_PcG_scaling.py $1 $2
