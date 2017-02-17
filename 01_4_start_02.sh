#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd
#$ -N linking_parts
#$ -l h_rt=0:00:10
#$ -l h_vmem=1G
#$ -pe sharedmem 1

. /etc/profile.d/modules.sh

module load anaconda/2.3.0
source activate Hi-C

SF=$(pwd)
ssh headnode1.ecdf.ed.ac.uk '
cd '$SF'
python 02_0_parse_merge_all.py' 
