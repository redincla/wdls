#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=Phasing
#SBATCH -t 4-00:00
#SBATCH --mem-per-cpu=10G
#SBATCH -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/logs/Phasing-sr.slurm.%N.%j.log

# checks for appropriate input parameters

source /dcsrsoft/spack/bin/setup_dcsrsoft
module load gcc
module load samtools/1.12

samtools reheader -c 'sed -e '28d' ' -P /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/Phasing/WES_Data/31691/test.bam > /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/Phasing/WES_Data/31691/test_clean.bam
#samtools view -H /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/Phasing/WES_Data/31691/test.bam | sed -e '28d' | samtools reheader - > /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/Phasing/WES_Data/31691/test_clean.bam
samtools index /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/Phasing/WES_Data/31691/test_clean.bam