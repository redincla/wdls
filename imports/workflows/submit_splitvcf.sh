#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=SplitVCF
#SBATCH -t 4-00:00
#SBATCH --mem-per-cpu=10G
#SBATCH -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/logs/splitVCF.slurm.%N.%j.log

# checks for appropriate input parameters
if [ $# -eq 2 ]; then
 input_file=$1 #full path to json input file for FilterVCF.wdl workflow
 options_file=$2 #full path to json options file for FilterVCF.wdl workflow

else
 echo -e "\n\nSubmit split and filter VCF workflow (FilterVCF.wdl) \n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "submit_split.vcf.sh [input_file] [option_file] "
 echo "input_file: full path to input json file"
 echo "option_file: full path to options json file"
 exit 1
fi

module add Development/java/1.8.0_232

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar run /home/credin/scratch/WGS/wdls/imports/workflows/FilterVCF.wdl -i ${input_file} -o ${options_file}
