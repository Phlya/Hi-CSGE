#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd
#$ -l h_rt=6:00:00
#$ -l h_vmem=64G
#$ -pe sharedmem 1
. /etc/profile.d/modules.sh

module load anaconda/2.3.0
module load igmm/libs/ncurses/6.0
module load igmm/apps/hdf5/1.8.16
source activate Hi-C

python 02_6_export_for_juicebox.py $1

