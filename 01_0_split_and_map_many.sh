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
    qsub 01_1_launch_splitting_and_mapping.sh $barcode $side 100000000
  done
done

sleep 1h
qsub -hold_jid splitting,mapping 01_4_start_02.sh
