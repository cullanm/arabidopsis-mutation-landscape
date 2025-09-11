#!/bin/bash
#SBATCH --job-name=duplex       			# Job name
#SBATCH --partition=schmitz_p               	# Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=1                  # Number of CPU cores per task
#SBATCH --mem=2G                          	# Job memory request snakemake rarely uses more than a ratio of 2G/core
#SBATCH --time=7-00:00:00                     # Time limit hrs:min:sec
#SBATCH --output=duplex.out         # Standard output log
#SBATCH --error=duplex.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cullan.meyer@uga.edu   	# Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node
SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2;exit}')  # location of this script
HOME_DIR=$(dirname $(dirname $SCRIPT_PATH))     # directory with the git repo (contains all the scripts)
SCRATCH_DIR=/scratch/cam02551/$(basename $HOME_DIR)     # directory with the data
set -euo pipefail

# load modules
ml Mamba

source activate snakemake	# this needs to be a mamba environment created from workflow/envs/snakemake.yaml

snakemake --profile ${HOME_DIR}/smk_profiles/slurm --directory ${SCRATCH_DIR} duplex_mutations duplex_stats
