#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd
#$ -l h_rt=6:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 8
. /etc/profile.d/modules.sh

module load java/jdk/1.8.0

cat /exports/eddie/scratch/s1529682/processed/merged/for_juicebox/$1\_for_juicebox.txt | sort -k2,2d -k6,6d -n -t $'\t' -T /exports/eddie/scratch/s1529682/tmp/ -S 128G --parallel 8 | uniq | nl -n ln > /exports/eddie/scratch/s1529682/processed/merged/for_juicebox/$1\_for_juicebox_sorted.txt
cd /exports/igmm/eddie/wendy-lab/ilia/scripts
sh juicebox.sh pre /exports/eddie/scratch/s1529682/processed/merged/for_juicebox/$1\_for_juicebox_sorted.txt $1.hic mm9
