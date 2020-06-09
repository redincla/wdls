#!/bin/bash

#SBATCH -n 1
#SBATCH --job-name=CollectMetrics
#SBATCH -t 0-03:00
#SBATCH --mem-per-cpu=14G
#SBATCH --output=slurm.%N.%j.log

module add Development/java/1.8.0_232
module add R/3.6.1

echo "hello from `hostname` at `date`"

java -Dconfig.file=/home/credin/.cromwell.conf -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/imports/tasks/CollectRawWgs-test.wdl -i /home/credin/scratch/WGS/wdls/imports/tasks/CollectRawWgs-test.json