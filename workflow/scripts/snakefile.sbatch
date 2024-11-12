#!/bin/bash
#SBATCH --job-name=smk_main
#SBATCH --output=smk_main_%j.out
#SBATCH --error=smk_main_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --partition=long
#SBATCH --mem=10G

#Unload any active module
module purge
#Load modules for snakemake (in the draco cluster)
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh
#Activate snakemake conda environment
conda activate snakemake_v7.24.0

#Run Snakemake
#snakemake --use-conda --conda-frontend conda --cores 1 --configfile config/config.json
snakemake --use-conda --conda-frontend conda --configfile config/config.json --jobs 50 --profile simple

#Deactivate conda environment
conda deactivate
