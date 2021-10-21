#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=extractVCF
#SBATCH -t 01-00:00
#SBATCH --mem-per-cpu=20G
#SBATCH -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/json/extractVCF.slurm.%N.%j.log

sample_ID=$1
module add Development/java/1.8.0_232
module add UHTS/Analysis/HTSlib/1.10.1
module add UHTS/Analysis/BEDTools/2.29.2
module add UHTS/Analysis/samtools/1.10

java -Xmx8g -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar \
    SelectVariants \
    -R /home/credin/refs/references/hg38/Homo_sapiens_assembly38.fasta \
    -V /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/210708_TACG_CoPRAC_HEV_VZVrefined_priors-filtered.vcf.gz \
    -O /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/Phasing/${sample_ID}.vcf.gz \
    -sn ${sample_ID} \
    --exclude-non-variants