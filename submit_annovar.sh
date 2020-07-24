#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=Annovar
#SBATCH -t 00-2:00
#SBATCH --mem-per-cpu=16G
#SBATCH -o slurm.%N.%j.log

module add Development/java/1.8.0_232
module add R/3.6.1

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/imports/workflows/annotateVCF.wdl -i /home/credin/scratch/WGS/wdls/imports/config-files/annotateVCF.inputs.json -o /home/credin/scratch/WGS/wdls/imports/config-files/annotateVCF.options.json