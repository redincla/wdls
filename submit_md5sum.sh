#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --job-name=md5sum
#SBATCH -t 720
#SBATCH --mem-per-cpu=1G
#SBATCH --output=slurm.%N.%j.log

#SBATCH --chdir /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/

echo "hello from `hostname` at `date`"
md5sum -c testCoPrAC_md5.txt 
