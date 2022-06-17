#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=Extractvcf
#SBATCH -t 00-12:00
#SBATCH --mem-per-cpu=10G
#SBATCH -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/logs/220302_extractVCFAli.slurm.%N.%j.log

module add Development/java/1.8.0_232

java -Xmx8g -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar \
SelectVariants \
-R /users/credin/refs/references/hg38/Homo_sapiens_assembly38.fasta \
-V /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/220204_278WGSrefined_priors-filtered.vcf.gz \
-O /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/tmp.vcf.gz \
--exclude-non-variants \
--restrict-alleles-to BIALLELIC \
-select "AF > 0.02"

java -Xmx8g -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar \
SelectVariants \
-R /users/credin/refs/references/hg38/Homo_sapiens_assembly38.fasta \
-V /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/tmp.vcf.gz \
-O /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/220204_278WGSrefined_priors-filtered-MAF2PCT.vcf.gz \
--exclude-non-variants \
--restrict-alleles-to BIALLELIC \
-select "AC > 4"

java -Xmx8g -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar \
VariantsToTable \
-V /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/220204_278WGSrefined_priors-filtered-MAF2PCT.vcf.gz \
-F CHROM -F POS -F REF -F ALT -F AC -F AN -F AF -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F VAR -F NSAMPLES -F NCALLED \
-O /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/220204_278WGSrefined_priors-filtered-MAF2PCT.tsv

    