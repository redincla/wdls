#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=SplitVCF
#SBATCH -t 00-04:00
#SBATCH --mem-per-cpu=10G
#SBATCH -o slurm.%N.%j.log

#sample_ID=$1
module add Development/java/1.8.0_232

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar run /home/credin/scratch/WGS/wdls/imports/workflows/FilterVCF.wdl -i /home/credin/scratch/WGS/wdls/imports/config-files/FilterVCF.inputs.json -o /home/credin/scratch/WGS/wdls/imports/config-files/FilterVCF.options.json