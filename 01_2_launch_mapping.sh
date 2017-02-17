#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd
#$ -N mapping
#$ -l h_rt=12:00:00
#$ -l h_vmem=5G
#$ -pe sharedmem 4

 sleep $[ ( $RANDOM % 60 )  + 1 ]s

. /etc/profile.d/modules.sh

module load anaconda/2.3.0
module load igmm/libs/ncurses/6.0
module load igmm/apps/bowtie/2.2.6
module load igmm/apps/samtools/1.3
module load igmm/apps/hdf5/1.8.16
source activate Hi-C

python 01_3_iterative_mapping.py $1
