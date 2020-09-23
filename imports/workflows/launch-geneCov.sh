 #!/bin/bash

module add Development/java/1.8.0_232
module add UHTS/Analysis/BEDTools/2.29.2

# checks for appropriate input
if [ $# -eq 3 ]; then
    sample_list=$1 #list of bam files to process
    bedfile=$2
	genome_index=$3
else
 echo -e "\n\nLaunching script to compute mean coverage by bedtools intervals by sample bam\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "launch-geneCov.sh [sample_list] [bedfile] [genome_index]"
 echo "sample_list: list of bam files to process"
 echo "bedfile: sorted bedfile with intervals to compute mean coverage on"
 echo "genome_index: reference genome index file"
 exit 1
fi

while read bam; do
sample_name=$(echo "${bam}" | cut -f 18 -d/ | cut -f 1 -d.)
sbatch -J geneCov_${sample_name} -D /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/json --time 240 -p normal \
  		-n 1 \
  		--mem-per-cpu=1G \
 	    --wrap "bedtools coverage -a ${bedfile} -b ${bam} -g ${genome_index} -sorted -mean -nobuf > /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/${sample_name}/QC/${sample_name}_cardioGenes_codingExons_coverage.list"
done < ${sample_list}