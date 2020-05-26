#!/bin/bash

#SBATCH -n 1
#SBATCH --job-name=WGS
#SBATCH -t 3-00:00
#SBATCH --mem-per-cpu=14G
#SBATCH --output=slurm.%N.%j.log
#module load Utility/cromwell/47

module add Development/java/1.8.0_232
module add R/3.6.1
module add UHTS/Analysis/samtools/1.10
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module add UHTS/Analysis/picard-tools/2.21.8

echo "hello from `hostname` at `date`"

java -Dconfig.file=/home/credin/.cromwell.conf -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/WholeGenomeGermlineSingleSample.wdl -i /home/credin/scratch/WGS/wdls/WholeGenomeGermlineSingleSample.json
