#!/bin/bash
#SBATCH --job-name=scrambler       			# Job name
#SBATCH --partition=schmitz_p               	# Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=8                  # Number of CPU cores per task
#SBATCH --mem=16G                          	# Job memory request snakemake rarely uses more than a ratio of 2G/core
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec
#SBATCH --output=scrambler_%a.out         # Standard output log
#SBATCH --error=scrambler_%a.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cullan.meyer@uga.edu   	# Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node
#SBATCH --array=6-6
SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2;exit}')  # location of this script
HOME_DIR=$(dirname $(dirname $SCRIPT_PATH))     # directory with the git repo (contains all the scripts)
SCRATCH_DIR=/scratch/cam02551/$(basename $HOME_DIR)     # directory with the data
cd ${SCRATCH_DIR}
set -euo pipefail

ml matplotlib/3.5.2-foss-2021b
ml SAMtools/1.14-GCC-11.2.0
ml Pysam/0.17.0-GCC-11.2.0
ml tqdm/4.62.3-GCCcore-11.2.0
ml SciPy-bundle/2021.10-foss-2021b

python3.9 ${HOME_DIR}/python_scripts/duplex_strand_swapper.py -@ 8 --replicate_chars 10 -o data/align/big_Col-0-swapped-${SLURM_ARRAY_TASK_ID}.bam data/align/big_Col-0-*_merged.bam
samtools index -@ 8 data/align/big_Col-0-swapped-${SLURM_ARRAY_TASK_ID}.bam

#python3.9 python_scripts/duplex_scrambled_subsampler.py -@ 8 --replicate_chars 10 -o data/align/big_Col-0-scrambled-1 data/align/big_Col-0-*_merged.bam
