#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=GatherSampleEvidence
#SBATCH -t 5-00:00
#SBATCH --mem-per-cpu=10G
#SBATCH -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/logs/GatherSampleEvidence.slurm.%N.%j.log

source /dcsrsoft/spack/bin/setup_dcsrsoft
module load  openjdk/16.0.1

if [ $# -eq 2 ]; then
 input_file=$1 #full path to json input file for runSmoove-largecohort.wdl workflow
 options_file=$2 #full path to json options file for runSmoove-largecohort.wdl workflow

else
 echo -e "\n\nSubmit GatherSampleEvidence workflow (GatherSampleEvidence.wdl) \n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "sbatch launch_GatherSampleEvidence.sh [input_file] [options_file] "
 echo "input_file: full path to inputs json file"
 echo "options_file: full path to options json file"
 exit 1
fi

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar run /home/credin/scratch/WGS/wdls/imports/workflows/GatherSampleEvidence.wdl -i ${input_file} -o ${options_file} 
