#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=runDelly
#SBATCH -t 07-00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/logs/slurm.Delly-FG.%N.%j.log

module add Development/java/1.8.0_232

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/imports/workflows/runDelly-sample.wdl -i /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/json/runDelly.inputs.json -o /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/json/runDelly.options.json