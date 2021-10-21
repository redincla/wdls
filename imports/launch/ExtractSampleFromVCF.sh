#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=VCFSplit
#SBATCH -t 00-08:00
#SBATCH --mem-per-cpu=9G
#SBATCH -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/VZV/Analysis/json/SplitBySample.slurm.%N.%j.log

sample_ID=$1
input_vcf=$2

module add Development/java/1.8.0_232

java -Xmx9g -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar \
     SelectVariants \
        -R /home/credin/refs/references/hg38/Homo_sapiens_assembly38.fasta \
        -V ${input_vcf} \
        -O "/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/${sample_ID}.vcf.gz" \
        --keep-original-ac \
        --exclude-non-variants \
        -sn ${sample_ID} -sn "28" -sn "29"

