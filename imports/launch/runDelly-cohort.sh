#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=runDelly
#SBATCH -t 07-00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/logs/slurm.Delly-FG.%N.%j.log

module add Development/java/1.8.0_232

if [ $# -eq 2 ]; then
 input_file=$1 #full path to json input file 
 options_file=$2 #full path to json options

else
 echo -e "\n\nSubmit Delly workflow on WGS cohort (runDelly-cohort.wdl) \n\n"
 echo "Usage:"
 echo "sbatch runDelly-cohort.sh [input_file] [options_file] "
 echo "input_file: full path to inputs json file"
 echo "options_file: full path to options json file"
 exit 1
fi


java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/imports/workflows/runDelly-cohort.wdl -i ${input_file} -o ${options_file} 