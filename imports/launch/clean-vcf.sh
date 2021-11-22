#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=clean-VCF
#SBATCH -t 01-00:00
#SBATCH --mem-per-cpu=30G
#SBATCH -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/logs/cleanVCF.slurm.%N.%j.log

module add UHTS/Analysis/EPACTS/3.2.6
module add UHTS/Analysis/samtools/1.10

cd /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointAnnotation

#bcftools view -h 211105_256WGSrefined_priors-filtered.annotated-FIX.vcf.gz > header
#gzip -cd 211105_256WGSrefined_priors-filtered.annotated-FIX.vcf.gz | grep -vP '^#' | sed 's@ @-@g' > 211105_256WGSrefined_priors-filtered.annotated-FIX-NO-HEADER.vcf

cat header 211105_256WGSrefined_priors-filtered.annotated-FIX-NO-HEADER.vcf | bgzip -c > 211105_256WGSrefined_priors-filtered.annotated-FIX-CLEAN.vcf.gz
tabix 211105_256WGSrefined_priors-filtered.annotated-FIX-CLEAN.vcf.gz