#!/bin/bash

#SBATCH -n 1
#SBATCH --job-name=WGS
#SBATCH -t 3-00:00
#SBATCH --mem-per-cpu=14G
#SBATCH --output=slurm.%N.%j.log
#module load Utility/cromwell/47

module add Development/java/1.8.0_232
module add R/3.6.1

echo "hello from `hostname` at `date`"

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/WholeGenomeGermlineSingleSample.wdl -i /home/credin/scratch/WGS/wdls/WholeGenomeGermlineSingleSample.json
