#!/bin/bash

full_list=$1

while read sample_ID; do

echo "#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=encrypt
#SBATCH -t 2-00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/json/${sample_ID}.encrypt.slurm.%N.%j.log

#cd /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/${sample_ID}/Raw
cd /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs
#cat ${sample_ID}*R1_001.fastq.gz > /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/${sample_ID}_R1_001.fq.gz
#cat ${sample_ID}*R2_001.fastq.gz > /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/${sample_ID}_R2_001.fq.gz
md5sum /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/${sample_ID}_R1_001.fq.gz >> /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/batch6_md5sum_file
md5sum /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/${sample_ID}_R2_001.fq.gz >> /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/batch6_md5sum_file
openssl aes-256-cbc -e -in /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/${sample_ID}_R1_001.fq.gz -out /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/crypted/${sample_ID}_R1_001_crypted.fq.gz -pass file:/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/keyfile_batch_6
openssl aes-256-cbc -e -in /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/${sample_ID}_R2_001.fq.gz -out /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/crypted/${sample_ID}_R2_001_crypted.fq.gz -pass file:/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/RFU_FQs/keyfile_batch_6" > /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/json/${sample_ID}.encrypt.script-submit

sbatch /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/json/${sample_ID}.encrypt.script-submit


done < ${full_list}
