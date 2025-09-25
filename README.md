# Description
Processing pipeline for nanorate sequencing data (NanoSeq) in [Nanorate sequencing reveals the _Arabidopsis_ somatic mutation landscape](https://doi.org/10.1101/2025.06.15.659769). This pipeline is intended to be usable on any Duplex sequencing data and run on either a SLURM computing cluster or locally.

# Running the pipeline

### 1. Create a directories for the code and data
```
mkdir {code dir}
mkdir {data dir}
```

### 2. Clone the repository
```
cd {code dir}
git clone https://github.com/cullanm/arabidopsis-mutation-landscape.git\ .
```

### 3. Install Mamba
Either by following the installation instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or loading it as a module on the cluster.

### 4. Create a Mamba environment for Snakemake
```
mamba env create -f workflow/envs/snakemake.yaml
```

### 5. Modify and run snakemake.sh
Open the script `sub_scripts/snakemake.sh` and set the HOME_DIR and SCRATCH_DIR variables to {code dir} and {data dir} respectively. Comment out SCRIPT_PATH. If running locally, change the `--profile` parameter of the snakemake command to `${HOME_DIR}/smk_profiles/local`

If running locally:
```
chmod +x sub_scripts/snakemake.sh
sub_scripts/snakemake.sh
```

If running on a SLURM cluster:

`sbatch sub_scripts/snakemake.sh`

# Running with your own data
To run the pipeline with different Duplex-seq data, do the following before step 5 above:

### Alter `sample_table.tsv`
Delete or comment out the existing rows of the table (apart from the header). Add a new row for each sequencing library. If using data on SRA, enter an SRA accession in the table and it will be downloaded as part of the pipeline. If using local data, move the fastq files into `{data dir}/data/fastq/` and rename them to `{group}_{subgroup}_{rep}.fastq{.gz if gzipped}` where {group} {subgroup} and {rep} are defined in `sample_table.tsv`. Then follow step 5 above.

# Output
This will produce the following output:

`{data dir}/data/variant/{group}_{subgroup}_merged_filtered.tsv`: filtered set of variants. This still includes variants with abundance >35% in the sample, which are likely germline. To get only somatic mutations, these should be filtered out using the sample_sup and sample_cov columns. The meaning of each column can be found in the `--help` messages of `python_scripts/duplex_caller.py` and `python_scripts/add_duplex_filter_columns.py`.

`{data dir}/data/coverage/{group}_{subgroup}_merged_duplex/Chr{1-5}.npy`: numpy arrays containing the callable coverage of every site in the genome (e.g. if the 10th position of the `Chr2.npy` array has a value of 5, that means there were 5 opportunities to detect a mutation at Chr2:10 (zero based) in that sample)

`{data dir}/data/metadata/duplex_metadata/metadata.tsv`: table with information on the number of reads and molecules in the sample (used for Table S4)

`{data dir}/data/metadata/duplex_metadata/conflict_rate.svg`: estimated conflict rates and imputed conflict rates (Figure S6)

`{data dir}/data/metadata/duplex_metadata/reads_per_molecule.svg`: histogram of reads per molecule (Figure S7)

`{data dir}/data/metadata/duplex_metadata/subsampling_efficiency.svg`: sequencing efficiencies (callable molecules / read pair) of libraries when subsampled (Figure S8)

`{data dir}/data/metadata/duplex_filters.svg`: visualization of how estimated false positive and true positive rates change as each filter threshold is adjusted (Figure S3)

`{data dir}/logs/`: contains printed output of all commands run. Some of these contain useful information (e.g. `filter_duplex_variants_{group}_{subgroup}_{rep}.err` contains the information of Table S5)

# Producing figures
Most of the figures were generated in `lab_notebook/mutation_landscape_figures.ipynb`. Before running this notebook, you'll need to run the whole pipeline and run/follow the instructions in `lab_notebook/reformatting_ma_line_muts.ipynb`. If you're unable to get the pipeline running, you can download the unfiltered variants and callable coverage at this [Zenodo repository](). The unfiltered variants must first be filtered using `python_scripts/filter_duplex_variants.py` with default parameters.

# Command list
If Snakemake isn't working or you don't want to use it. Here's a list of the commands needed to produce the output. Only commands for sample untreated_1_a are shown here. These commands should all be run from {data dir}.
```
# download the reference genome
wget --no-verbose -O data/ref/ref.fa.gz https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz &> logs/fetch_ref_fasta.out
gunzip data/ref/ref.fa.gz

# make a bowtie2 index
mkdir -p data/ref/bowtie_index &> null || true
bowtie2-build --threads 1 data/ref/ref.fa data/ref/bowtie_index/index &> logs/ref_bowtie_index.out

# fetch from sra
prefetch SRR33007850 --max-size 100GB -O data/fastq &> logs/fetch_srr_SRR33007850.out
fasterq-dump data/fastq/SRR33007850 -O data/fastq -t tmp/ -e 4 &> logs/extract_srr_paired_SRR33007850.out
mv data/fastq/SRR33007850_1.fastq data/fastq/untreated_1_a_1.fastq

# trim reads
fastp -i data/fastq/untreated_1_a_1.fastq -I data/fastq/untreated_1_a_2.fastq -o data/fastq/untreated_1_a_1_trim.fastq -O data/fastq/untreated_1_a_2_trim.fastq -h logs/trim_reads_paired_untreated_1_a.html -w 16  &> logs/trim_reads_paired_untreated_1_a.out

# align reads
bowtie2 -x data/ref/bowtie_index/index -1 data/fastq/untreated_1_a_1_trim.fastq -2 data/fastq/untreated_1_a_2_trim.fastq -S data/align/untreated_1_a_bowtie.sam -p 16 --time -X 800 &> logs/align_fastq_bowtie_untreated_1_a.out

# sort, mark duplicates, and filter reads (note: filtering of optical duplicates will not occur if using data from SRA, as read names are overwritten. This may slightly alter the output. Reach out to the corresponding author if original reads are needed)
samtools fixmate -@ 8 -u -m data/align/untreated_1_a_bowtie.sam data/align/untreated_1_a_fixmate.bam &> logs/sam_fixmate_untreated_1_a.out
sambamba sort -t 12 -m 512MiB --compression-level 1 -o data/align/untreated_1_a_sort.bam data/align/untreated_1_a_fixmate.bam &> logs/sam_sort_untreated_1_a.out
samtools markdup -@ 8  -d 2500 -s data/align/untreated_1_a_sort.bam data/align/untreated_1_a_markdup.bam &> logs/sam_markdup_untreated_1_a.out
sambamba view -t 8 --format bam -l 3 --filter "mapping_quality >= 1 and proper_pair and ([dt] == null or [dt] != 'SQ')" -o data/align/untreated_1_a_filtered.bam data/align/untreated_1_a_markdup.bam

# merge replicate libraries (above code must have finished for all samples before running)
samtools merge -@ 8 -r -o data/align/untreated_1_merged.bam data/align/untreated_1_a_filtered.bam data/align/untreated_1_b_filtered.bam data/align/untreated_1_c_filtered.bam &> logs/merge_replicate_bams_untreated_1.out
samtools index -@ 8 data/align/untreated_1_merged.bam &> logs/index_bam_untreated_1_merged.out

# calculate genome coverage
bedtools genomecov -bg -ibam data/align/untreated_1_merged.bam > data/coverage/untreated_1_merged.bg 2> logs/coverage_bg_untreated_1_merged.out

# make duplex blacklists (above code must have finished for all samples before running)
python {code dir}/python_scripts/generate_duplex_blacklists.py -@ 8 --fasta data/ref/ref.fa --bedgraphs data/coverage/untreated_1_merged.bg data/coverage/untreated_2_merged.bg data/coverage/untreated_3_merged.bg data/coverage/untreated_4_merged.bg data/coverage/untreated_5_merged.bg data/coverage/untreated_6_merged.bg data/coverage/untreated_7_merged.bg data/coverage/untreated_8_merged.bg --output data/ref/duplex_blacklist/ &> logs/generate_duplex_blacklists.out

# make variant tables
python {code dir}/python_scripts/duplex_caller.py -@ 8 --umi RG data/align/untreated_1_merged.bam data/variant/untreated_1_merged_duplex.tsv &> logs/duplex_caller_untreated_1_merged.out
python {code dir}/python_scripts/dup_informed_caller.py -@ 8 --umi RG --duplicate_support 2 data/align/untreated_1_merged.bam data/variant/untreated_1_merged_informed.tsv &> logs/dup_informed_caller_untreated_1_merged.tsv

# filter variants (above code must have finished for all samples before running). Note that untreated_1 is not used as a control for itself
python {code dir}/python_scripts/add_duplex_filter_columns.py --input data/variant/untreated_1_merged_duplex.tsv --frequency data/variant/untreated_1_merged_informed.tsv --controls data/variant/untreated_2_merged_informed.tsv data/variant/untreated_3_merged_informed.tsv data/variant/untreated_4_merged_informed.tsv data/variant/untreated_5_merged_informed.tsv data/variant/untreated_6_merged_informed.tsv data/variant/untreated_7_merged_informed.tsv data/variant/untreated_8_merged_informed.tsv --blacklists data/ref/duplex_blacklist/ --output data/variant/untreated_1_merged_added.tsv &> logs/add_duplex_filter_columns_untreated_1_merged.out
python {code dir}/python_scripts/filter_duplex_variants.py --input data/variant/untreated_1_merged_added.tsv --output data/variant/untreated_1_merged_filtered.tsv > logs/filter_duplex_variants_untreated_1_merged.out 2> logs/filter_duplex_variants_untreated_1_merged.err

# calculate callable coverage
python {code dir}/python_scripts/duplex_coverage.py -@ 8 --umi RG --blacklist data/ref/duplex_blacklist/ data/align/untreated_1_merged.bam data/coverage/untreated_1_merged_duplex/ &> logs/duplex_coverage_untreated_1_merged.out

# make metadata table and figures
python {code dir}/python_scripts/duplex_metadata.py -@ 8 --umi RG --replicate ,0 --region Chr1:0-1000000 --input data/align/untreated_1_merged.bam data/align/untreated_2_merged.bam data/align/untreated_3_merged.bam data/align/untreated_4_merged.bam data/align/untreated_5_merged.bam data/align/untreated_6_merged.bam data/align/untreated_7_merged.bam data/align/untreated_8_merged.bam --output data/metadata/duplex_metadata/ > logs/duplex_metadata.out 2> logs/duplex_metadata.err

# make swapped samples, call variants, and calculate callable coverage (repeat each of these commands 4 times, making duplex_swapped_{0-3})
python {code dir}/python_scripts/duplex_strand_swapper.py -@ 8 --umi RG --input data/align/untreated_1_merged.bam data/align/untreated_2_merged.bam data/align/untreated_3_merged.bam data/align/untreated_4_merged.bam data/align/untreated_5_merged.bam data/align/untreated_6_merged.bam data/align/untreated_7_merged.bam data/align/untreated_8_merged.bam --output data/align/duplex_swapped_0.bam &> logs/duplex_strand_swapper_0.out
samtools index -@ 1 data/align/duplex_swapped_0.bam &> logs/index_bam_duplex_swapped_0.out
python {code dir}/python_scripts/duplex_caller.py -@ 8 --umi RG data/align/duplex_swapped_0.bam data/variant/duplex_swapped_0_duplex.tsv &> logs/duplex_caller_duplex_swapped_0.out
python {code dir}/python_scripts/dup_informed_caller.py -@ 8 --umi RG --duplicate_support 2 data/align/duplex_swapped_0.bam data/variant/duplex_swapped_0_informed.tsv &> logs/dup_informed_caller_duplex_swapped_0.tsv
python {code dir}/python_scripts/add_duplex_filter_columns.py --input data/variant/duplex_swapped_0_duplex.tsv --frequency data/variant/duplex_swapped_0_informed.tsv --controls data/variant/untreated_1_merged_informed.tsv data/variant/untreated_2_merged_informed.tsv data/variant/untreated_3_merged_informed.tsv data/variant/untreated_4_merged_informed.tsv data/variant/untreated_5_merged_informed.tsv data/variant/untreated_6_merged_informed.tsv data/variant/untreated_7_merged_informed.tsv data/variant/untreated_8_merged_informed.tsv --blacklists data/ref/duplex_blacklist/ --output data/variant/duplex_swapped_0_added.tsv &> logs/add_duplex_filter_columns_duplex_swapped_0.out
python {code dir}/python_scripts/duplex_coverage.py -@ 8 --umi RG --blacklist data/ref/duplex_blacklist/ data/align/duplex_swapped_0.bam data/coverage/duplex_swapped_0_duplex/ &> logs/duplex_coverage_duplex_swapped_0.out


# generate filter threshold false positive rate estimates
python {code dir}/python_scripts/visualize_duplex_filters.py --real_vars data/variant/untreated_1_merged_added.tsv data/variant/untreated_2_merged_added.tsv data/variant/untreated_3_merged_added.tsv data/variant/untreated_4_merged_added.tsv data/variant/untreated_5_merged_added.tsv data/variant/untreated_6_merged_added.tsv data/variant/untreated_7_merged_added.tsv data/variant/untreated_8_merged_added.tsv --swapped_vars data/variant/duplex_swapped_0_added.tsv data/variant/duplex_swapped_1_added.tsv data/variant/duplex_swapped_2_added.tsv data/variant/duplex_swapped_3_added.tsv --real_cov data/coverage/untreated_1_merged_duplex data/coverage/untreated_2_merged_duplex data/coverage/untreated_3_merged_duplex data/coverage/untreated_4_merged_duplex data/coverage/untreated_5_merged_duplex data/coverage/untreated_6_merged_duplex data/coverage/untreated_7_merged_duplex data/coverage/untreated_8_merged_duplex/ --swapped_cov data/coverage/duplex_swapped_0_duplex data/coverage/duplex_swapped_1_duplex data/coverage/duplex_swapped_2_duplex data/coverage/duplex_swapped_3_duplex/ --output data/metadata/duplex_filters.svg &> logs/visualize_duplex_filters.out
```
