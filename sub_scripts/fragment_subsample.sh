#!/bin/bash
#SBATCH --job-name=fragment_subsample       			# Job name
#SBATCH --partition=schmitz_p               	# Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=8                  # Number of CPU cores per task
#SBATCH --mem=16G                          	# Job memory request snakemake rarely uses more than a ratio of 2G/core
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec
#SBATCH --output=fragment_subsample_%a.out         # Standard output log
#SBATCH --error=fragment_subsample_%a.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cullan.meyer@uga.edu   	# Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node
#SBATCH --array=1-8
SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2;exit}')  # location of this script
HOME_DIR=$(dirname $(dirname $SCRIPT_PATH))     # directory with the git repo (contains all the scripts)
SCRATCH_DIR=/scratch/cam02551/$(basename $HOME_DIR)     # directory with the data
cd ${SCRATCH_DIR}
set -euo pipefail

ml SciPy-bundle/2022.05-foss-2022a
ml tqdm/4.64.0-GCCcore-11.3.0
ml Pysam/0.19.1-GCC-11.3.0
ml matplotlib/3.5.2-foss-2022a
ml SAMtools/1.16.1-GCC-11.3.0

########################################## code for subsampling the NanoSeq libraries to match the sequencing depth of deduplicated ones
#cd data/align
#mkdir -p wgs_matched_subsample
#file=big_Col-0-${SLURM_ARRAY_TASK_ID}_merged.bam
#declare -A sub_fracs=( [1]=0.38579 [2]=0.22570 [3]=0.23019 [4]=0.24623 [5]=0.35991 [6]=0.13354 [7]=0.24452 [8]=0.27696)
#
#python3.10 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 --umi RG --random_frac ${sub_fracs[$SLURM_ARRAY_TASK_ID]} -o wgs_matched_subsample/$(basename -s _merged.bam $file)_matched.bam $file
#samtools index wgs_matched_subsample/$(basename -s _merged.bam $file)_matched.bam

########################################## subsample NanoSeq libraries to measure mutation rate stability
mkdir -p data/align/fragment_subsample
cd data/align/fragment_subsample
file=big_Col-0-${SLURM_ARRAY_TASK_ID}_merged.bam
python3.10 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 --umi RG --random_frac 0.8 -o $(basename -s _merged.bam $file)_08frags.bam ../$file
python3.10 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 --umi RG --random_frac 0.6 -o $(basename -s _merged.bam $file)_06frags.bam ../$file
python3.10 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 --umi RG --random_frac 0.4 -o $(basename -s _merged.bam $file)_04frags.bam ../$file
python3.10 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 --umi RG --random_frac 0.2 -o $(basename -s _merged.bam $file)_02frags.bam ../$file
python3.10 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 --umi RG --random_frac 0.1 -o $(basename -s _merged.bam $file)_01frags.bam ../$file

samtools index -@ 6 $(basename -s _merged.bam $file)_08frags.bam
samtools index -@ 6 $(basename -s _merged.bam $file)_06frags.bam
samtools index -@ 6 $(basename -s _merged.bam $file)_04frags.bam
samtools index -@ 6 $(basename -s _merged.bam $file)_02frags.bam
samtools index -@ 6 $(basename -s _merged.bam $file)_01frags.bam

########################################## code for subsampling NanoSeq libraries to test if low coverage fragments are prone to false positives
#cd data/align
#files=(big_Col-0-*_merged.bam)
#file=${files[SLURM_ARRAY_TASK_ID]}
#base=$(basename -s _merged.bam $file)

#mkdir -p fragment_subsample
#python3.9 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 --required_strand_coverage 8 --umi RG -o fragment_subsample/${base}_8s_16t_nocap.bam $file
#samtools index fragment_subsample/${base}_8s_16t_nocap.bam

#python3.9 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 -s 6       --umi RG -o fragment_subsample/${base}_6s_12t.bam ${base}_8s_16t_nocap.bam
#python3.9 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 -s 4 -t 10 --umi RG -o fragment_subsample/${base}_4s_10t.bam ${base}_8s_16t_nocap.bam
#python3.9 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 -s 4       --umi RG -o fragment_subsample/${base}_4s_8t.bam ${base}_8s_16t_nocap.bam
#python3.9 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 -s 2 -t 6  --umi RG -o fragment_subsample/${base}_2s_6t.bam ${base}_8s_16t_nocap.bam
#python3.9 ${HOME_DIR}/python_scripts/duplex_subsampler.py -@ 6 -s 2       --umi RG -o fragment_subsample/${base}_2s_4t.bam ${base}_8s_16t_nocap.bam
#
#samtools index -@ 6 fragment_subsample/${base}_6s_12t.bam
#samtools index -@ 6 fragment_subsample/${base}_4s_10t.bam
#samtools index -@ 6 fragment_subsample/${base}_4s_8t.bam
#samtools index -@ 6 fragment_subsample/${base}_2s_6t.bam
#samtools index -@ 6 fragment_subsample/${base}_2s_4t.bam
