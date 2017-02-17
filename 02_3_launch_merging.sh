#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd
#$ -l h_rt=6:00:00
#$ -l h_vmem=32G
#$ -pe sharedmem 4
. /etc/profile.d/modules.sh

module load anaconda/2.3.0
module load igmm/libs/ncurses/6.0
module load igmm/apps/hdf5/1.8.16
source activate Hi-C

LC_ALL=C
python 02_4_merge_fragments.py $1

ssh headnode1.ecdf.ed.ac.uk "cd /exports/igmm/eddie/wendy-lab/ilia/scripts; qsub -N export"$1" 02_5_launch_export.sh "$1""
ssh headnode1.ecdf.ed.ac.uk "cd /exports/igmm/eddie/wendy-lab/ilia/scripts; qsub -N sorting"$1" -hold_jid export"$1" 02_7_sort_for_juicebox.sh "$1""

for b in 1000 100 25 10 5
do
  cooler cload --hiclib --assembly mm9 $b\Kb_bins_mm9.bed /exports/eddie/scratch/s1529682/processed/merged/$1\_fragment_dataset.hdf5 /exports/eddie/scratch/s1529682/processed/merged/$1\_$b\Kb.cool
  cooler balance -p 4 /exports/eddie/scratch/s1529682/processed/merged/$1\_$b\Kb.cool
done
