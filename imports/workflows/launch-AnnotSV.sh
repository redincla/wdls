#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=annotSV
#SBATCH -t 2-00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=AnnotSV.slurm.%N.%j.log
    
export ANNOTSV=/home/credin/refs/tools/AnnotSV  #add ANNOTSV to path
module add UHTS/Analysis/samtools/1.10  #load bcftools
module load UHTS/Analysis/BEDTools/2.29.2 #load bedtools

$ANNOTSV/bin/AnnotSV -SvinputFile /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/Delly/Full_Genome/TACG_COPRAC_HEV_VZV_0621.delly.vcf.gz \
    -genomeBuild GRCh38 \
    -minTotalNumber 100 \
    -outputFile /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/Delly/Full_Genome/210701_AnnotSV/TACG_COPRAC_HEV_VZV_0621.delly.annotSV_full+PASS.tsv \
    -candidateGenesFile /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/Delly/Full_Genome/CAR_IFN.list \
    -candidateGenesFiltering no \
#    -candidateSnvIndelFiles /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/210617_TACG_CoPRAC_HEV_VZVrefined_priors-filtered.vcf.gz \
    -snvIndelFiles /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/210617_TACG_CoPRAC_HEV_VZVrefined_priors-filtered.vcf.gz \
    -snvIndelPASS 1 \
#    -rankFiltering "3,4,5" \    