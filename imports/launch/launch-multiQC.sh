#!/bin/bash

# checks for appropriate input
if [ $# -eq 2 ]; then
    input_dir=$1 #path to analysis directory 
    output_dir=$2 #path to output directory for html reports
    echo "#!/bin/bash
    #SBATCH -n 1
    #SBATCH --job-name=multiQC
    #SBATCH -t 2-00:00
    #SBATCH --mem-per-cpu=1G
    #SBATCH --output=${output_dir}/multiQC.slurm.%N.%j.log
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc python
    source /home/credin/refs/tools/multiQC-venv/bin/activate
    multiqc ${input_dir} --outdir ${output_dir}" > ${output_dir}/multiQC.script-submit

elif [ $# -eq 3 ]; then
    file_list=$2 #optional, list of QC files to process, if spreadout/not all within the same folder
    output_dir=$3 #path to output directory for html reports
    echo "#!/bin/bash
    #SBATCH -n 1
    #SBATCH --job-name=multiQC
    #SBATCH -t 2-00:00
    #SBATCH --mem-per-cpu=1G
    #SBATCH --output=${output_dir}/multiQC.slurm.%N.%j.log
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc python
    source /home/credin/refs/tools/multiQC-venv/bin/activate
    multiqc --file-list ${file_list} --outdir ${output_dir}" > ${output_dir}/multiQC.script-submit
   
else
 echo -e "\n\nLaunching multiQC\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "launch-multiQC.sh [input_dir] [output_dir] (default)"
 echo "launch-multiQC.sh F [file_list] [output_dir]"
 echo "input_dir: directory to parse QC files"
 echo "file_list: alternative argument, TO BE PRECEDED WITH 'F', file with list of all QC reports to parse, when they are spread out"
 echo "output_dir: path where to drop multiQC reports"
 exit 1
fi

sbatch ${output_dir}/multiQC.script-submit