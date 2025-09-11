#!/bin/bash
#SBATCH --job-name=read_subsample       			# Job name
#SBATCH --partition=schmitz_p               	# Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=16                  # Number of CPU cores per task
#SBATCH --mem=24G                          	# Job memory request snakemake rarely uses more than a ratio of 2G/core
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec
#SBATCH --output=read_subsample_%a.out         # Standard output log
#SBATCH --error=read_subsample_%a.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cullan.meyer@uga.edu   	# Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node
#SBATCH --array=0-7
SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2;exit}')  # location of this script
HOME_DIR=$(dirname $(dirname $SCRIPT_PATH))     # directory with the git repo (contains all the scripts)
SCRATCH_DIR=/scratch/cam02551/$(basename $HOME_DIR)     # directory with the data
cd ${SCRATCH_DIR}
set -euo pipefail

ml Sambamba

cd data/align
files=(big_Col-0-*_merged.bam)
file=${files[SLURM_ARRAY_TASK_ID]}
base=$(basename -s .bam $file)

mkdir -p read_subsample
#sambamba view -t 16 -h -s 0.8 -f bam -l 1 -o read_subsample/${base}_08.bam $file
#sambamba view -t 16 -h -s 0.6 -f bam -l 1 -o read_subsample/${base}_06.bam $file
#sambamba view -t 16 -h -s 0.4 -f bam -l 1 -o read_subsample/${base}_04.bam $file
#sambamba view -t 16 -h -s 0.2 -f bam -l 1 -o read_subsample/${base}_02.bam $file
sambamba view -t 16 -h -s 0.1 -f bam -l 1 -o read_subsample/${base}_01.bam $file
