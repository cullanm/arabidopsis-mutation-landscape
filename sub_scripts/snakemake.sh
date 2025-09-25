#!/bin/bash
#SBATCH --job-name=snakemake       			# Job name
#SBATCH --partition=schmitz_p               	# Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=1                  # Number of CPU cores per task
#SBATCH --mem=2G                          	# Job memory request snakemake rarely uses more than a ratio of 2G/core
#SBATCH --time=7-00:00:00                     # Time limit hrs:min:sec
#SBATCH --output=snakemake.out         # Standard output log
#SBATCH --error=snakemake.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cullan.meyer@uga.edu   	# Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node
SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2;exit}')  # COMMENT OUT IF NOT RUNNING ON SLURM
HOME_DIR=$(dirname $(dirname $SCRIPT_PATH))     # SET TO THE DIRECTORY OF THE GIT REPO (CODE DIR)
SCRATCH_DIR=/scratch/cam02551/$(basename $HOME_DIR)     # SET TO THE DIRECTORY OF THE DATA (DATA DIR)
set -euo pipefail

# load modules
ml Mamba

source activate snakemake	# this needs to be a mamba environment created from workflow/envs/snakemake.yaml

snakemake --profile ${HOME_DIR}/smk_profiles/slurm --directory ${SCRATCH_DIR} --verbose duplex_mutations duplex_stats
