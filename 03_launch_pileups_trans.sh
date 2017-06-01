#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd
#$ -N pileups
#$ -l h_rt=6:00:00
#$ -l h_vmem=36G
#$ -pe sharedmem 4
. /etc/profile.d/modules.sh

module load anaconda/2.3.0
module load igmm/apps/hdf5/1.8.16
source activate Hi-C

python 03_pileups_trans.py $1 $2 40
