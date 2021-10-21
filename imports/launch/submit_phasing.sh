#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=Phasing
#SBATCH -t 4-00:00
#SBATCH --mem-per-cpu=10G
#SBATCH -o /scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/logs/Phasing-sr.slurm.%N.%j.log

# checks for appropriate input parameters
if [ $# -eq 2 ]; then
 input_file=$1 #full path to json input file for Phasing-sr.wdl workflow
 options_file=$2 #full path to json options file for Phasing-sr.wdl workflow

else
 echo -e "\n\nSubmit phasing for Illumina short reads workflow (Phasing-sr.wdl) \n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "submit_phasing.sh [input_file] [option_file] "
 echo "input_file: full path to input json file"
 echo "option_file: full path to options json file"
 exit 1
fi

module add Development/java/1.8.0_232

java -Dconfig.file=/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/data_and_refs/data/Raw_FQ/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar run /scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/workflows/Phasing-sr.wdl -i ${input_file} -o ${options_file}
