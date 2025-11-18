#!/bin/bash
#SBATCH --job-name=snakemake		# Job name
#SBATCH --partition=batch		# Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1			# Run a single task
#SBATCH --cpus-per-task=1		# Number of CPU cores per task (increase if running with the local profile)
#SBATCH --mem=2G			# Job memory request (increase to ~2G/core if running with the local profile)
#SBATCH --time=7-00:00:00		# Time limit hrs:min:sec
#SBATCH --output=snakemake.out		# Standard output log
#SBATCH --error=snakemake.err		# Standard error log
#SBATCH --export=NONE			# do not load any env variables to compute node
HOME_DIR=$(dirname $(dirname $SCRIPT_PATH))     # SET TO THE DIRECTORY OF THE GIT REPO (CODE DIR)
SCRATCH_DIR=/scratch/cam02551/$(basename $HOME_DIR)     # SET TO THE DIRECTORY OF THE DATA (DATA DIR)
set -euo pipefail

# load modules (if running on cluster)
ml Mamba

source activate snakemake	# this needs to be a mamba environment created from workflow/envs/snakemake.yaml

# for running on a slurm cluster
snakemake --profile ${HOME_DIR}/smk_profiles/slurm --directory ${SCRATCH_DIR} --verbose duplex_mutations duplex_stats

# for running locally
#snakemake --profile ${HOME_DIR}/smk_profiles/local --directory ${SCRATCH_DIR} --cores 64 --resources mem_mb=131072 --verbose duplex_mutations duplex_stats
