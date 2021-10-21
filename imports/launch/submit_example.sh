#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --job-name=TO-FILL
#SBATCH -t 0-24:00
#SBATCH --mem-per-cpu=1G
#SBATCH --output=slurm.%N.%j.log

#SBATCH --chdir ~/.

echo "hello from `hostname` at `date`"
