#!/bin/bash

module add Development/java/1.8.0_232
module add UHTS/Quality_control/fastqc/0.11.7


# checks for appropriate input
if [ $# -eq 1 ]; then
    sample_list=$1 #list of bam files to process
else
 echo -e "\n\nLaunching fastQC script by sample bam\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "launch-fastQC.sh [sample_list] "
 echo "sample_list: list of bam files to process"
 exit 1
fi

while read bam; do
sbatch -J fastQC -D /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/json --time 1400 -p normal \
  		-n 1 \
  		--mem-per-cpu=8G \
 	    --wrap "fastqc ${bam} -o /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/fastQC -f bam"
done < ${sample_list}