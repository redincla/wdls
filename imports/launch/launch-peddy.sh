#!/bin/bash

####################################
#    launch Peddy    #
####################################

# Copyright (c) 2020 Claire Redin and The CHUV Precision Medicine Unit, Lausanne, Switzerland.
# Contact: Claire Redin <claire.redin@chuv.ch>
# Peddy was developped by the Quinlan lab: https://github.com/brentp/peddy 

# checks for appropriate input
if [ $# -eq 3 ]; then
 vcf_file=$1  
 ped_file=$2
 OUTDIR=$3

else
 echo -e "\n\nLaunching peddy to perform relatedness and sex check from WGS data \n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "launch-peddy.sh [gvcf] [ped_file] [OUTDIR] "
 echo "vcf_file: full path to multi-sample gvcf file"
 echo "ped_file: full path to cohort ped file"
 echo "OUTDIR: full path to output directory"
 exit 1
fi

# Load Peddy virtual environment
source /dcsrsoft/spack/bin/setup_dcsrsoft
module load gcc python
source /home/credin/refs/tools/peddy_venv/bin/activate

cd ${OUTDIR}

### Write and launch slurm script
sbatch -J Peddy -D ${OUTDIR} --time 200 -p normal \
  		-n 1 \
  		--mem-per-cpu=2G \
 	    --wrap "peddy --plot -p 3 --sites hg38 --prefix HEV ${vcf_file} ${ped_file}"