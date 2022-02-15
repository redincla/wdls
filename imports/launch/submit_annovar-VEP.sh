#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=JointAnnotation
#SBATCH -t 3-00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/logs/slurm.JointAnnotation-VEP.%N.%j.log

module add Development/java/1.8.0_232

if [ $# -eq 2 ]; then
 input_file=$1 #full path to json input file for annotateVCF.wdl workflow
 options_file=$2 #full path to json options file for annotateVCF.wdl workflow

else
 echo -e "\n\nSubmit annotation workflow (annotateVCF.wdl) \n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "sbatch submit_annovar-VEP.sh [input_file] [options_file] "
 echo "input_file: full path to inputs json file"
 echo "options_file: full path to options json file"
 exit 1
fi

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar run /home/credin/scratch/WGS/wdls/imports/workflows/annotateVCF-VEP.wdl -i ${input_file} -o ${options_file} 
