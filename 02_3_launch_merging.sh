#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 8
. /etc/profile.d/modules.sh

module load anaconda/2.3.0
module load igmm/libs/ncurses/6.0
module load igmm/apps/hdf5/1.8.16
source activate Hi-C

LC_ALL=C
python 02_4_merge_fragments.py $1

for b in 1000 100 25 10 5
do
  cooler cload --hiclib --assembly mm9 $b\Kb_bins_mm9.bed /exports/eddie/scratch/s1529682/processed/merged/$1\_fragment_dataset.hdf5 /exports/eddie/scratch/s1529682/processed/merged/$1\_$b\Kb.cool
  cooler balance -p 8 /exports/eddie/scratch/s1529682/processed/merged/$1\_$b\Kb.cool
done

python 02_5_export_for_juicebox.py $1
cat /exports/eddie/scratch/s1529682/processed/merged/for_juicebox/$1\_for_juicebox.txt | sort -k2,2d -k6,6d -n -t $'\t' -T /exports/eddie/scratch/s1529682/tmp/ -S 128G --parallel 8 | uniq | nl -n ln > /exports/eddie/scratch/s1529682/processed/merged/for_juicebox/$1\_for_juicebox_sorted.txt
