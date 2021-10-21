#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=testVEP
#SBATCH -t 4-00:00
#SBATCH --mem-per-cpu=10G
#SBATCH -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/logs/VEP.slurm.%N.%j.log

input_vcf=$1

source /dcsrsoft/spack/bin/setup_dcsrsoft
module load singularity
export PERL5LIB=$PERL5LIB:/db/local/vep/Plugins.  # point to plugin repository
export SINGULARITY_BINDPATH="/scratch,/users,/dcsrsoft,/db"

singularity run /dcsrsoft/singularity/containers/ensembl-vep_104.sif vep \
    -i ${input_vcf} \
    --plugin SpliceAI --dir_plugins /db/local/vep/plugins_data/ \
    --buffer_size 100000 \
    --offline --fork 10 \
    --dir_cache=/db/local \
    --vcf \
    --cache -o /users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/JointCalling/Trio_27-28-29.VEPspliceAI.vcf.gz