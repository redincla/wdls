 #!/bin/bash

module add Development/java/1.8.0_232
module add UHTS/Analysis/BEDTools/2.29.2

# checks for appropriate input
if [ $# -eq 2 ]; then
    sample_list=$1 #list of bam files to process
    bedfile=$2
else
 echo -e "\n\nLaunching script to compute mean coverage by bedtools intervals by sample bam\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "launch-geneCov.sh [sample_list] [bedfile]"
 echo "sample_list: list of bam files to process"
 echo "bedfile: sorted bedfile with intervals to compute mean coverage on"
 exit 1
fi

while read bam; do
sbatch -J geneCov -D /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/json --time 120 -p normal \
  		-n 1 \
  		--mem-per-cpu=1G \
 	    --wrap "bedtools coverage -a ${bedfile} -b ${bam} -sorted -mean -nobuf > test
done < ${sample_list}
