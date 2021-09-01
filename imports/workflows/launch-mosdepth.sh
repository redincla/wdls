#!/bin/bash

full_list=$1

while read sample_ID cram; do
echo "#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=mosdepth
#SBATCH -t 2-00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/json/${sample_ID}.mosdepth.slurm.%N.%j.log

source /dcsrsoft/spack/bin/setup_dcsrsoft
module load gcc python
source /software/etc/profile.d/conda.sh
conda activate /home/credin/refs/tools/.conda/envs/mosdepth_env
mosdepth --fast-mode -f /home/credin/refs/references/hg38/Homo_sapiens_assembly38.fasta ${sample_ID} ${cram} -b /home/credin/refs/references/gene_lists/cardioGenes_full-CDS+coding_exons_faisorted.bed -T 5,10,20 -n" > /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/json/${sample_ID}.mosdepth.script-submit

cd /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/${sample_ID}/QC/
sbatch /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/json/${sample_ID}.mosdepth.script-submit

done < ${full_list}