#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=annotSV
#SBATCH -t 2-00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --output=AnnotSV.slurm.%N.%j.log
    
export ANNOTSV=/home/credin/refs/tools/AnnotSV  #add ANNOTSV to path
module add UHTS/Analysis/samtools/1.10  #load bcftools
module load UHTS/Analysis/BEDTools/2.29.2 #load bedtools

input_vcf=$1
output_file=$2

$ANNOTSV/bin/AnnotSV -SvinputFile ${input_vcf} \
    -genomeBuild GRCh38 \
    -minTotalNumber 100 \
    -outputFile ${output_file} \
    -candidateGenesFile /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/Delly/Full_Genome/CAR_IFN.list \
    -candidateGenesFiltering no \
    -snvIndelFiles /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/220204_278WGSrefined_priors-filtered.vcf.gz \
    -snvIndelPASS 1 \
    -rankFiltering "3,4,5"
    #    -candidateSnvIndelFiles /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/210617_TACG_CoPRAC_HEV_VZVrefined_priors-filtered.vcf.gz \
