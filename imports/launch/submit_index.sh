#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=index_dict
#SBATCH -t 0-02:00
#SBATCH --mem-per-cpu=6G
#SBATCH --output=slurm.%N.%j.log

module add Development/java/1.8.0_232
cd /users/credin/refs/references/GRCh38
java -jar /software/UHTS/Analysis/picard-tools/2.21.8/bin/picard.jar CreateSequenceDictionary R=GCA_000001405.15_GRCh38_no_alt_analysis_set.fa O=GCA_000001405.15_GRCh38_no_alt_analysis_set.dict