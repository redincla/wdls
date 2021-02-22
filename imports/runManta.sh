#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=runManta
#SBATCH -t 04-00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/logs/slurm.Manta-FG.%N.%j.log

module add Development/java/1.8.0_232

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/imports/workflows/runManta.wdl -i /home/credin/scratch/WGS/wdls/imports/config-files/Manta.inputs.json -o /home/credin/scratch/WGS/wdls/imports/config-files/Manta.options.json