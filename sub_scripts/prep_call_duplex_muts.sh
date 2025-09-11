#!/bin/bash
#SBATCH --job-name=prep_call_duplex_muts       			# Job name
#SBATCH --partition=batch               	# Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=8                          # Run a single task	
#SBATCH --cpus-per-task=6                  # Number of CPU cores per task
#SBATCH --mem-per-cpu=2G                          	# Job memory request snakemake rarely uses more than a ratio of 2G/core
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec
#SBATCH --output=prep_call_duplex_muts.out         # Standard output log
#SBATCH --error=prep_call_duplex_muts.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cullan.meyer@uga.edu   	# Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node
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
# also need genmap

shopt -s extglob
files=(data/align/big_Col-0-+([0-9])_merged.bam)
bases=()
for f in ${files[@]}; do
        echo $f
        bases+=($(basename -s .bam $f))
done
echo ${bases[@]}

# count coverage
# notes on parellelizing with srun: --export=ALL must be set to inherit loaded modules, --nodes must be 1 or srun will try to take more for some reason, --exclusive must be set or some sruns won't start, --mem-per-cpu must be specified instead of --mem or the first srun might take all the memory
#srun --ntasks=1 --nodes=1 --export=ALL --exclusive bedtools genomecov -bg -ibam ${files[0]} > tmp/${bases[0]}.bg &
#srun --ntasks=1 --nodes=1 --export=ALL --exclusive bedtools genomecov -bg -ibam ${files[1]} > tmp/${bases[1]}.bg &
#srun --ntasks=1 --nodes=1 --export=ALL --exclusive bedtools genomecov -bg -ibam ${files[2]} > tmp/${bases[2]}.bg &
#srun --ntasks=1 --nodes=1 --export=ALL --exclusive bedtools genomecov -bg -ibam ${files[3]} > tmp/${bases[3]}.bg &
#srun --ntasks=1 --nodes=1 --export=ALL --exclusive bedtools genomecov -bg -ibam ${files[4]} > tmp/${bases[4]}.bg &
#srun --ntasks=1 --nodes=1 --export=ALL --exclusive bedtools genomecov -bg -ibam ${files[5]} > tmp/${bases[5]}.bg &
#srun --ntasks=1 --nodes=1 --export=ALL --exclusive bedtools genomecov -bg -ibam ${files[6]} > tmp/${bases[6]}.bg &
#srun --ntasks=1 --nodes=1 --export=ALL --exclusive bedtools genomecov -bg -ibam ${files[7]} > tmp/${bases[7]}.bg &
#wait
#echo done counting coverage


# make control files using dup_informed_caller
srun --ntasks=1 --nodes=1 --export=ALL --exclusive python ${HOME_DIR}/python_scripts/dup_informed_caller.py --umi RG --duplicate_support 2 -@ 6 ${files[0]} data/variant/${bases[0]}_informed.tsv &
srun --ntasks=1 --nodes=1 --export=ALL --exclusive python ${HOME_DIR}/python_scripts/dup_informed_caller.py --umi RG --duplicate_support 2 -@ 6 ${files[1]} data/variant/${bases[1]}_informed.tsv &> null &
srun --ntasks=1 --nodes=1 --export=ALL --exclusive python ${HOME_DIR}/python_scripts/dup_informed_caller.py --umi RG --duplicate_support 2 -@ 6 ${files[2]} data/variant/${bases[2]}_informed.tsv &> null &
srun --ntasks=1 --nodes=1 --export=ALL --exclusive python ${HOME_DIR}/python_scripts/dup_informed_caller.py --umi RG --duplicate_support 2 -@ 6 ${files[3]} data/variant/${bases[3]}_informed.tsv &> null &
srun --ntasks=1 --nodes=1 --export=ALL --exclusive python ${HOME_DIR}/python_scripts/dup_informed_caller.py --umi RG --duplicate_support 2 -@ 6 ${files[4]} data/variant/${bases[4]}_informed.tsv &> null &
srun --ntasks=1 --nodes=1 --export=ALL --exclusive python ${HOME_DIR}/python_scripts/dup_informed_caller.py --umi RG --duplicate_support 2 -@ 6 ${files[5]} data/variant/${bases[5]}_informed.tsv &> null &
srun --ntasks=1 --nodes=1 --export=ALL --exclusive python ${HOME_DIR}/python_scripts/dup_informed_caller.py --umi RG --duplicate_support 2 -@ 6 ${files[6]} data/variant/${bases[6]}_informed.tsv &> null &
srun --ntasks=1 --nodes=1 --export=ALL --exclusive python ${HOME_DIR}/python_scripts/dup_informed_caller.py --umi RG --duplicate_support 2 -@ 6 ${files[7]} data/variant/${bases[7]}_informed.tsv &> null &
wait

# make blacklist arrays
python ${HOME_DIR}/python_scripts/generate_duplex_blacklists.py -@ 8 --fasta data/ref/ref.fa --output data/blacklist/arabidopsis/ --bedgraphs \
tmp/${bases[0]}.bg \
tmp/${bases[1]}.bg \
tmp/${bases[2]}.bg \
tmp/${bases[3]}.bg \
tmp/${bases[4]}.bg \
tmp/${bases[5]}.bg \
tmp/${bases[6]}.bg \
tmp/${bases[7]}.bg

echo completed prep_call_duplex_muts.sh
