#!/bin/bash

#SBATCH -n 1
#SBATCH --job-name=aggregateBamQC
#SBATCH -t 0-12:00
#SBATCH --mem-per-cpu=10G

#SBATCH -o slurm.%N.%j.log
#module load Utility/cromwell/47


module add Development/java/1.8.0_232
module add R/3.6.1

echo "hello from `hostname` at `date`"

java -Dconfig.file=/home/credin/.cromwell.conf -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/imports/workflows/AggregatedBamQC.wdl -i /home/credin/scratch/WGS/wdls/imports/workflows/AggregatedBamQC.json