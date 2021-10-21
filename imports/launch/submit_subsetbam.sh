#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=SubsetBam
#SBATCH -t 00-04:00
#SBATCH --mem-per-cpu=9G
#SBATCH -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/json/subsetCram.slurm.%N.%j.log

sample_ID=$1
module add Development/java/1.8.0_232
module add UHTS/Analysis/HTSlib/1.10.1
module add UHTS/Analysis/BEDTools/2.29.2
module add UHTS/Analysis/samtools/1.10

export CRAM_REFERENCE="/home/credin/refs/references/hg38/Homo_sapiens_assembly38.fasta"

#bedtools intersect -a /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/${sample_ID}/Processed/${sample_ID}.cram \
#-b /home/credin/refs/references/gene_lists/cardioGenes_txn_canonical_annotated.NO-HEADER_faisorted.bed -sorted > /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/subset_cram/${sample_ID}_CardioGenes.cram
#samtools index /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/subset_cram/${sample_ID}_CardioGenes.cram

#bedtools intersect -a /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/${sample_ID}/Processed/${sample_ID}.cram \
#-b /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/subset_bam/HLA_alleles_faisorted.bed > /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/subset_bam/${sample_ID}_HLA.cram
#samtools index /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/subset_bam/${sample_ID}_HLA.cram

samtools view /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/Phasing/${sample_ID}.cram "chr11" -C \
    -T /home/credin/refs/references/hg38/Homo_sapiens_assembly38.fasta > /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/Phasing/${sample_ID}.chr11.cram
samtools index /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/CoPrAc/2101/Phasing/${sample_ID}.chr11.cram

