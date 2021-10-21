#!/bin/bash

echo "#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=DoC
#SBATCH -t 05:00
#SBATCH --mem-per-cpu=12G
#SBATCH --output=/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/json/DoC.slurm.%N.%j.log

module add Development/java/1.8.0_232
module add UHTS/Analysis/BEDTools/2.29.2

bedtools intersect -a ~/refs/references/GRCh38/exome_coverage/all.90pctOver10X.LiftedTohg38.bed   \
-b /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L1/QC/L1.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_08/QC/LUG_08.per-base.bed.gz \
-c > /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/BurdenTest/MAF1%/test.DoC.bed" > /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/json/intersect.script-submit
#/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L1/QC/L1.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L10/QC/L10.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L11/QC/L11.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L12/QC/L12.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L13/QC/L13.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L2/QC/L2.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L3/QC/L3.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L4/QC/L4.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L5/QC/L5.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L6/QC/L6.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L7/QC/L7.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L8/QC/L8.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/L9/QC/L9.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_02/QC/LUG_02.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_04/QC/LUG_04.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_05/QC/LUG_05.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_06/QC/LUG_06.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_07/QC/LUG_07.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_08/QC/LUG_08.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_09/QC/LUG_09.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_10/QC/LUG_10.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_11/QC/LUG_11.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_12/QC/LUG_12.per-base.bed.gz /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/LUG_13/QC/LUG_13.per-base.bed.gz

sbatch /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/HEV/Analysis/json/intersect.script-submit