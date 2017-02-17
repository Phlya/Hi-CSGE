mkdir -p /exports/eddie/scratch/s1529682/fastq/to_process
mkdir -p /exports/eddie/scratch/s1529682/fastq/split
mkdir -p /exports/eddie/scratch/s1529682/bams
mkdir -p /exports/eddie/scratch/s1529682/tmp/
mkdir -p /exports/eddie/scratch/s1529682/processed/merged/for_juicebox/

for barcode in "$@"
do
  echo $barcode
  for side in 1 2
  do
    echo $side
    qsub 01_1_launch_splitting_and_mapping.sh $barcode $side 40000000
    qsub 01_1_launch_splitting_and_mapping.sh $barcode\t $side 40000000 #For test sequencing files with different read length, so can't just combine them
  done
done

#sleep 1m
#qsub -hold_jid splitting,mapping 01_4_start_02.sh
