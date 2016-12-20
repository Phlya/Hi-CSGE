#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd
#$ -N splitting
#$ -l h_rt=03:00:00
#$ -l h_vmem=2G
#$ -pe sharedmem 4
. /etc/profile.d/modules.sh
module load igmm/apps/pigz/2.3.3
SF=$(pwd)
cd /exports/eddie/scratch/s1529682/fastq/to_process

pigz -dc -p 4 $1_$2.fq.gz | split -l $3 - -d -a 3 "../split/$1_$2.fq.gz"
ssh headnode1.ecdf.ed.ac.uk '
cd '$SF'
for file in /exports/eddie/scratch/s1529682/fastq/split/'$1_$2'.fq.gz*
do
  echo $file
  qsub 01_2_launch_mapping.sh $file
done'

