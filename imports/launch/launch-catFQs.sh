#!/bin/bash

full_list=$1

while read sample_ID; do

echo "#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=cat_fq
#SBATCH -t 2-00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/json/${sample_ID}.cat-FQ.slurm.%N.%j.log

#cd /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/${sample_ID}/Raw/
cd /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/WGS_Fellay_TACG_Apr21
cat ${sample_ID}*R1_001.fastq.gz > /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/${sample_ID}_R1_001.fq.gz
cat ${sample_ID}*R2_001.fastq.gz > /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/${sample_ID}_R2_001.fq.gz" > /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/json/${sample_ID}.cat-FQ.script-submit

sbatch /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/json/${sample_ID}.cat-FQ.script-submit

done < ${full_list}