#!/bin/bash
#SBATCH --job-name=call_duplex_muts       			# Job name
#SBATCH --partition=batch               	# Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=6                  # Number of CPU cores per task
#SBATCH --mem=32G                          	# Job memory request snakemake rarely uses more than a ratio of 2G/core
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec
#SBATCH --output=call_duplex_muts_%a.out         # Standard output log
#SBATCH --error=call_duplex_muts_%a.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cullan.meyer@uga.edu   	# Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node
#SBATCH --array=0-7
SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2;exit}')  # location of this script
HOME_DIR=$(dirname $(dirname $SCRIPT_PATH))     # directory with the git repo (contains all the scripts)
SCRATCH_DIR=/scratch/cam02551/$(basename $HOME_DIR)     # directory with the data
cd ${SCRATCH_DIR}
set -euo pipefail

# load modules
ml matplotlib/3.9.2-gfbf-2024a
ml tqdm/4.66.5-GCCcore-13.3.0
ml Pysam/0.22.1-GCC-13.3.0
ml BEDTools/2.31.1-GCC-13.3.0
ml SAMtools/1.21-GCC-13.3.0

shopt -s extglob
#files=(data/align/big_@(arp|atx12r7|chr8|clf|ddm1|h2aw7|h2ax35|ku80|met1|msh2|msh6|nrpd1b|parp2|poly|rad5a|rad7a|sdg8|suvh456|xpg)*_merged.bam)
files=(data/align/big_Col-0-+([0-9])_merged.bam)
file=${files[SLURM_ARRAY_TASK_ID]}
base=$(basename -s .bam $file)

#samtools index -@ 6 $file

echo sample is $base

controls=""
for ((i=1;i<=8;i++)); do
	if [ $base != "big_Col-0-${i}_merged" ]; then
		controls+=" data/variant/big_Col-0-${i}_merged_informed.tsv"
	fi
done

echo $controls

# run duplex_caller, informed_caller for calculating variant frequencies, and informed_caller for control filters
#python ${HOME_DIR}/python_scripts/duplex_caller.py -@ 6 --umi RG $file data/variant/script_${base}_duplex.tsv 
#python ${HOME_DIR}/python_scripts/dup_informed_caller.py --umi RG --duplicate_support 2 -@ 6 $file data/variant/script_${base}_informed.tsv # these are used for determining the frequency of a mutation in a sample
#python ${HOME_DIR}/python_scripts/add_duplex_filter_columns.py --input data/variant/script_${base}_duplex.tsv --frequency data/variant/${base}_informed.tsv --controls $controls \
#--blacklists data/blacklist/arabidopsis/duplex_blacklist_ --output data/variant/script_${base}_added.tsv
python ${HOME_DIR}/python_scripts/filter_duplex_variants.py --input data/variant/script_${base}_added.tsv --output data/variant/script_${base}_filtered.tsv

# do a preliminary filter, require 1f1r original strand coverage and 24% support in each original strand
# duplex caller output is chrom pos ref alt frag_start frag_len frag_umi f_read1_sup f_read1_cov r_read1_sup r_read1_cov f_sup f_cov r_sup r_cov f_read1_conc r_read1_conc mq bq
#awk -F"\t" '{if ($1 == "chrom" || ($9 >= 1 && $11 >= 1 && $8 / $9 >= 0.24 && $10 / $11 >= 0.24)) print $0}' data/variant/${base}_duplex.tsv > data/variant/${base}_duplex_strandsup.tsv
#rm data/variant/${base}_duplex.tsv

# calculate duplex coverage
#mkdir -p data/coverage/2f1r
#python ${HOME_DIR}/python_scripts/duplex_coverage.py -@ 6 --umi RG -b data/blacklist/arabidopsis/duplex_blacklist_ $file data/coverage/2f1r/${base}_
