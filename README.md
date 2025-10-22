# Description
Processing pipeline for nanorate sequencing data (NanoSeq) in [Nanorate sequencing reveals the _Arabidopsis_ somatic mutation landscape](https://doi.org/10.1101/2025.06.15.659769). This pipeline can be run on either a SLURM computing cluster or locally and should work for any kind of Duplex Sequencing data (NanoSeq, BotSeqS, etc.) with or without in-read UMIs (though I recommend using UMIs if you're making libraries). All custom Duplex-seq python scripts can be found in the `python_scripts` directory and can be run directly from command line (details [below](#NanoSeq-data-processing-scripts)).

# Running the pipeline

## 1. Create a directories for the code and data
```
mkdir {code dir}
mkdir {data dir}
```

## 2. Clone the repository
```
cd {code dir}
git clone https://github.com/cullanm/arabidopsis-mutation-landscape.git .
```

## 3. Install Mamba
Either by following the installation instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or loading it as a module on your cluster.

## 4. Create a Mamba environment for Snakemake
```
mamba env create -f workflow/envs/snakemake.yaml
```

If this command fails, you can try setting your channel priority to flexible with `conda config --set channel_priority flexible` (or to 'strict' if your ~/.condarc was already set to flexible) or making the environment using `workflow/envs/snakemake_export.yaml` instead. If you still aren't able to make the environment work, you can try installing the packages within `workflow/envs/snakemake.yaml` using another method such as pip or loading cluster modules (note that if loading modules, the `module load` commands need to be present in both `sub_scripts/snakemake.sh` and `smk_profiles/slurm/slurm-jobscript.sh`).

## 5. Modify snakemake.sh (and settings.json)
Open the script `sub_scripts/snakemake.sh` and set the HOME_DIR and SCRATCH_DIR variables to {code dir} and {data dir} respectively. Comment out SCRIPT_PATH.

#### Running on a SLURM cluster
You may need to change the partition name in the header of `sub_scripts/snakemake.sh` and in `smk_profiles/slurm/settings.json`.

#### Running locally
If you're running the pipeline on a local machine or the cluster instructions did not work, use the alternative snakemake command in `sub_scripts/snakemake.sh`, specifying the number of cores and memory for snakemake to use in the command (and sbatch header).

## 6. Run Snakemake
#### Running on a SLURM cluster
```
sbatch sub_scripts/snakemake.sh
```

#### Running locally
```
chmod +x sub_scripts/snakemake.sh
sub_scripts/snakemake.sh
```

# Running with your own data
To run the pipeline with different Duplex-seq data, do the following before step 6 above:

### Alter `sample_table.tsv`
Delete or comment out the existing rows of the table (apart from the header). Add a new row for each sequencing library. If using data on SRA, enter an SRA accession in the table and it will be downloaded as part of the pipeline. If using local data, move the fastq files into `{data dir}/data/fastq/` and rename them to `{group}_{subgroup}_{rep}.fastq{.gz if gzipped}` where {group} {subgroup} and {rep} are defined in `sample_table.tsv`.

### Alter the config file
The file `config_files/arabidopsis.yaml` contains several settings for altering the pipeline's behavior, including specifying a different reference genome, different UMI structure (if not using xGen CS adapters), and different mutation filtering thresholds. If you're unsure of what filtering thesholds to use, I recommend running the pipeline with default parameters and viewing `{data dir}/data/metadata/duplex_filters.svg`. You can then make any changes to the filter thresholds in `config_files/arabidopsis.yaml`, delete all the filtered varaint files (`{data dir}/data/variant/*_filtered.tsv`), and rerun snakemake to regenerate them using the new thresholds.

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

# NanoSeq data processing scripts
Many of these scripts require numpy, pandas, tqdm, and pysam to be installed. `generate_duplex_blacklists.py` also requires genmap.

- [duplex_caller.py](#duplex_callerpy): makes an unfiltered table of variants, calulating the number of reads which support/cover variants on a molecule-by-molecule basis
- [dup_informed_caller.py](#dup_informed_callerpy): makes an unfiltered table of variants, calculating the number of molecules (not reads) which support/cover a variant
- [add_duplex_filter_columns.py](#add_duplex_filter_columnspy): adds additional information needed for filtering to the output of `duplex_caller.py`
- [filter_duplex_variants.py](#filter_duplex_variantspy): filters the output of `add_duplex_filter_columns.py` to a final set of mutations
- [duplex_coverage.py](#duplex_coveragepy): calculates the number of opportunities to detect a mutation at each site in the genome. Important for calculating mutation rates
- [duplex_strand_swapper.py](#duplex_strand_swapperpy): makes swapped null libraries out of duplex-seq libraries which can be used to estimate false positive rates
- [duplex_metadata.py](#duplex_metadatapy): reports statistics on the sequencing efficiency of duplex-seq libraries
- [visualize_duplex_filters.py](#visualize_duplex_filterspy): visualizes how adjusting the filters would alter true and false positive rates

### duplex_caller.py
```
usage: duplex_caller.py [-h] [-@ THREADS_INT] [-t TMP_STRING] [-u TAG] [-a] [-C INT] INPUT_FILE OUTPUT_FILE

Calls variants from a coordinate sorted BAM file of duplex sequencing reads. Ignores reads that didn't map as a concordant pair. Output will treat back-to-back SNVs and indels longer than 1bp as multiple 1bp variants, these will need
to be grouped later. For insertions, this means a CG insertion will have two entries in the output, a *->C and *->CG, with the second entry listing the support for the inserted G only. Note, the vast majority of variants output by
this script will not be true variants, further filtering steps will be necessary. Output is a tsv with the following columns:
- chrom: chromosome/contig
- pos: 0-based position of the variant. For insertions, this is the site in the reference which the insertion occurs before
- ref: reference sequence, "*" for insertion
- alt: alt sequence, "*" for deletion
- frag_start: start position of the original fragment which contains the variant (two variants at the same position may be on different fragments)
- frag_length: length of the original fragment which contains the variant
- frag_umi: UMI of the original fragment (will be empty if umi argument isn't specified)
- f_read1_sup: number of reads generated from the forward strand of the original fragment which support the variant
- f_read1_cov: number of reads generated from the forward strand of the original fragment which overlap the variant
- r_read1_sup: number of reads generated from the reverse strand of the original fragment which support the variant
- r_read1_cov: number of reads generated from the reverse strand of the original fragment which overlap the variant
- f_sup: number of reads mapping to the forward strand which support the variant
- f_cov: number of reads mapping to the forward strand which overlap the variant
- r_sup: number of reads mapping to the reverse strand which support the variant
- r_cov: number of reads mapping to the reverse strand which overlap the variant
- f_read1_conc: number of read pairs from the forward strand of the original fragment where both read1 and read2 support the variant
- r_read1_conc: number of read pairs from the reverse strand of the original fragment where both read1 and read2 support the variant
- mq: mapping quality of the reads supporting the variant. ascii converted (qual + 33)
- base quality: base call quality of the variant bases. ascii converted (qual + 33)
- end_mismatch_rate: maximum mismatch rate of any window starting from a fragment end and containing the variant. e.g. Removing all variants with end_mismatch_rate > 0.5 would mean ignoring the ends of each fragment until there are
more matches than mismatches from the position to the fragment end.

positional arguments:
  INPUT_FILE            BAM file to call from. Must be coordinate sorted and have an index file. Reads must be paired end with no unpaired reads present.
  OUTPUT_FILE           output variant file

options:
  -h, --help            show this help message and exit
  -@ THREADS_INT, --threads THREADS_INT
                        chromosomes are divided between threads. Cannot use more threads than num chromosomes + 1 (default: 1)
  -t TMP_STRING, --tmp TMP_STRING
                        prefix to append temporary output files (default: ./)
  -u TAG, --umi TAG     Consider the UMI in the provided BAM tag when grouping pre-PCR fragments (default: no UMI)
  -a, --all             when set, outputs all variants rather than only those with support in both original strands
  -C INT, --chunk_size INT
                        No effect on output. Number of reads to process before flushing buffered information. Increasing may speed up runtime, especially for high-coverage samples. Decreasing lowers memory requirement (default: 10000)
```

### dup_informed_caller.py
```
usage: dup_informed_caller.py [-h] [-@ INT] [-t PREFIX] [-q INT] [-d FLOAT] [-p INT] [-m INT] [-s INT] [-u TAG] [-e] [-C INT] FILE FILE

Calls variants from a coordinate sorted BAM file without any fancy realignment stuff. Ignores reads that didn't map as a concordant pair. Output will treat back-to-back SNVs as multiple 1bp variants. Indels of length >1 bp will appear
as a single row. Takes pcr/optical duplicates and read1-read2 overlap into consideration when determining if a read is in support and doesn't double count duplicates. For each set of reads with the same fragment start and end, they
only contribute to support if variants are supported by > x% of the covering reads in the set (includes both forward and reverse reads). If a set of reads passes this filter, it contributes 1 to forward/reverse support (or both if in
the read1-read2 overlap). This means read pairs with no duplicates still must pass this filter within the read1-read2 overlap. Outputs a tsv with columns:
- chrom: chromosome/contig
- pos: 0-based position of the variant. For insertions, this is the site in the reference which the insertion occurs before
- ref: reference sequence, "*" for insertion
- alt: alt sequence, "*" for deletion
- f_sup: number of forward reads supporting the variant
- f_cov: number of forward reads covering the variant position
- r_sup: number of reverse reads supporting the variant
- r_cov: number of reverse reads covering the variant position
- mapping quality: ascii converted (qual + 33) phred scaled mapping qualities of supporting reads. Length equals alt depth. Selects the median value for sets of duplicates
- base quality: ascii converted (qual + 33) phred scaled base call qualities of bases in supporting reads. Length equals alt depth * length of variant (e.g. a variant *->TA with two supporting reads where the base call qualities for
the two reads are ",," and "<<" will have a base quality of ",,<<"). Selects the median value for sets of duplicates

positional arguments:
  FILE                  BAM file to call from. Must be coordinate sorted and have an index file
  FILE                  output variant file

options:
  -h, --help            show this help message and exit
  -@ INT, --threads INT
                        chromosomes are divided between threads. Cannot use more threads than num chromosomes + 1 (default: 1)
  -t PREFIX, --tmp PREFIX
                        prefix to append temporary output files (default: ./)
  -q INT, --min_mq INT  minimum read mapping quality required for a read to be processed. Reads with lower mapping quality will not count towards any columns (default: 1)
  -d FLOAT, --duplicate_ratio FLOAT
                        Cutoff for whether a group of duplicates supports a variant. If the number of supporting reads from all duplicates / the number of covering reads from all duplicates (both forward and reverse) is greater than
                        DUPLICATE_RATIO, the group of duplicates will add 1 to the variant's support (1 to forward and reverse support if in the read1-read2 overlap) (default: 0.79)
  -p INT, --duplicate_support INT
                        Minimum number of duplicates required to contibute to variant support and coverage (default: 1)
  -m INT, --min_support INT
                        only output variants with f_sup + r_sup > MIN_SUPPORT (default: 1)
  -s INT, --min_strand_support INT
                        only output variants with f_sup > MIN_STRAND_SUPPORT and r_sup > MIN_STRAND_SUPPORT (default: 0)
  -u TAG, --umi TAG     Consider the UMI in the provided BAM tag when grouping pre-PCR fragments (default: no UMI)
  -e, --end_dist        Report the distance from fragment end of each supporting read as an extra column (default: not reported)
  -C INT, --chunk_size INT
                        No effect on output. Number of reads to process before flushing buffered information. Increasing may speed up runtime, especially for high-coverage samples. Decreasing lowers memory requirement (default: 10000)
```

### generate_duplex_blacklists.py
```
usage: generate_duplex_blacklists.py [-h] -f FASTA -b BAMs [BAMs ...] [-o PREFIX] [-l INT] [-m INT] [-@ INT] [-t PREFIX]

Creates blacklist files for filtering variants of duplex-seq libraries. Takes the reference genome and a set of BAM files for coverage calculations as input. Requires GenMap (https://github.com/cpockrandt/genmap) to be installed.
Produces four .npy arrays for each chromosome:
- {prefix}{chrom}_coverage_percentile.npy: the coverage percentile of each site in the genome, e.g. 0.01 means 1% of the genome has a lower or equal coverage
- {prefix}{chrom}_duplicate_distance.npy: the minimum hamming distance between a kmer at the genomic position of length --fragment_length and any other position in the genome. e.g. 2 means there's another kmer in the genome where only
2 bases differ. If no kmer is discovered, a value of 10 is assigned.
- {prefix}{chrom}_poly_at_length.npy: length of the longest poly-A/T repeat starting at the position or 1bp away. Only considers lengths >=4. e.g.
         sequence: ACTGAAAAAAAAGTCCCA
         blacklist: 000880000008800000
- {prefix}{chrom}_poly_repeat_length.npy: length of the longest poly-di/trinucleotide repeat starting at the position of 1bp away. Only considers length >=3. e.g.
         sequence: CCCGAGAGAGAGACCTC
         blacklist: 00554000000455000

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        genome fasta
  -b BAMs [BAMs ...], --bedgraphs BAMs [BAMs ...]
                        one or more coverage bedgraph files used to calculate coverage percentile per genomic site. These can be generated with "bedtools genomecov -bg -ibam {BAM}"
  -o PREFIX, --output PREFIX
                        prefix to prepend to blacklist outputs (default: ./)
  -l INT, --fragment_length INT
                        length of DNA fragments. Used for calculating the duplicate distance blacklist and passed to genmap -k argument (default: 150)
  -m INT, --max_hamming INT
                        maximum hamming distance to check with genmap for making the duplicate distance blacklist. Increasing this increases runtime. Setting to a negative value will skip running genmap and the blacklist will be all
                        10s (default: 3)
  -@ INT, --threads INT
                        number of threads to use when calculating genome coverage (default: 1)
  -t PREFIX, --tmp PREFIX
                        prefix to add to temporary files. Will make any directories which do not exist (default: tmp/)
```

### add_duplex_filter_columns.py
```
usage: add_duplex_filter_columns.py [-h] -i TSV -f TSV -c TSV [TSV ...] -b PREFIX [-o TSV]

Adds columns to the duplex-seq variants generated by duplex_caller.py which are required for filtering with filter_duplex_variants.py.
Requires the following as input:
- duplex variant tsv: generated by running duplex_caller.py on the sample to call variants within
- variant frequency tsv: generated by running dup_informed_caller.py --duplicate_support 2 on the same sample sample as the duplex variant tsv. This will be used to calculate the fraction of molecules in the sample which support each
variant
- control tsvs: generated by running dup_informed_caller.py --duplicate_support 2 --min_support 2 on different samples. This will be used to filter out any variants which occur too frequently in these samples
- blacklists: generated by generate_duplex_blacklists.py
Will add the following columns to the duplex variant tsv:
- high_sup: number of supporting bases with BQ>=30
- avg_mq: average MQ of supporting reads
- cov_per: coverage percentile of the variant site, uses blacklist (e.g. 0.02 means 2% of the genome has a lower read coverage than the variant position)
- dup_dist: minimum hamming distance between the fragment containing the variant and another genomic position, uses blacklist (e.g. 2 means there's another kmer in the genome that's 2 edits away from the kmer starting at the fragment
start (the value of "k" is defined when blacklists are generated and is not necessarily the fragment length))
- poly_at: length of longest poly-A/T repeat ending at or 1bp away from the variant position, uses blacklist (e.g. 5 means the variant is at the end of a 5bp poly-A or poly-T repeat)
- poly_rep: length of longest poly-di/trinucleotide repeat ending at or 1bp away from the variant position, uses blacklist (e.g. 4 means the variant is at the end of a length 3 dinucleotide or trinucleotide repeat (e.g. CCTCCTCCT))
- sample_sup: number of fragments in the sample which support the variant, uses the variant frequency tsv
- sample_cov: number of fragments in the sample which overlap the variant position, uses the variant frequency tsv
- control_sup: number of fragments in the control samples which support the variant, uses the control tsvs for this|n

options:
  -h, --help            show this help message and exit
  -i TSV, --input TSV   generated by running duplex_caller.py on the sample to call variants within. Must be sorted by chrom-pos-ref-alt, which duplex_caller.py does by default.
  -f TSV, --frequency TSV
                        generated by running dup_informed_caller.py --duplicate_support 2 on the same sample sample as --input. Must be sorted by chrom-pos-ref-alt, which dup_informed_caller.py does by default.
  -c TSV [TSV ...], --controls TSV [TSV ...]
                        generated by running dup_informed_caller.py --duplicate_support 2 --min_support 2 on samples other than --input
  -b PREFIX, --blacklists PREFIX
                        prefix of the blacklist files, same as was passed to generate_duplex_blacklists --prefix. Will look for {prefix}{chrom}_{blacklist type}.npy, where blacklist type is coverage_percentile, duplicate_distance,
                        poly_at_length, and poly_repeat_length.
  -o TSV, --output TSV  output tsv. Will have the same variants as --input but with extra columns added (default: stdout)
```

### filter_duplex_variants.py
```
usage: filter_duplex_variants.py [-h] [-i TSV] [-o TSV] [-A FLOAT] [-B INT] [-C INT] [-D INT] [-E FLOAT] [-F FLOAT] [-G INT] [-H FLOAT] [-I INT] [-J INT] [-K INT] [-L INT] [-M INT] [-N FLOAT] [-O INT] [-P FLOAT] [-Q]

Filters a duplex-seq variant tsv that has been generated by add_duplex_filter_columns.py. The output may still contain germline mutations and contaminant errors which will need to be filtered out by removing variants present in
multiple samples. The letter code of each filter is the same as in Figure S1. Filters are applied in alphabetical order.

options:
  -h, --help            show this help message and exit
  -i TSV, --input TSV   generated by running add_duplex_filter_columns.py on the output of duplex_caller.py (default: stdin)
  -o TSV, --output TSV  filtered output tsv. Will just be a subset of the input rows (default: stdout)
  -A FLOAT, --min_frac_init FLOAT
                        initial minimum fraction of supporting reads per strand (default: 0.24)
  -B INT, --min_strand_cov INT
                        minimum reads per strand overlapping variant (default: 2)
  -C INT, --min_total_cov INT
                        minimum total reads overalpping variant (default: 6)
  -D INT, --min_high_bq INT
                        minimum supporting reads with high BQ (>=30) (default: 4)
  -E FLOAT, --min_avg_mq FLOAT
                        minimum average MQ of supporting reads (default: 10)
  -F FLOAT, --max_end_mismatch FLOAT
                        maximum mismatch rate between a position and a fragment end before clipping (default: 0.21)
  -G INT, --min_indel_end_dist INT
                        minimum distance between indels and a fragment end (default: 6)
  -H FLOAT, --min_cov_percentile FLOAT
                        minimum sequencing coverage percentile of the variant position (default: 0.007, i.e. 0.7 percentile)
  -I INT, --max_poly_at INT
                        maximum length of a poly-A/T repeat ending at the variant position (default: 7)
  -J INT, --max_poly_rep INT
                        maximum length of a poly-di/trinucleotide repeat end at the variant position (default: 4)
  -K INT, --max_control_sup INT
                        maximum support for the variant across the control samples (default: 6)
  -L INT, --max_frag_len INT
                        maximum fragment length (default: 300)
  -M INT, --max_contamination_dist INT
                        maximum distance between two variants passing prior filters (A-L) within a single fragment (default: 4, i.e. there can be at most 3 bases between passing variants)
  -N FLOAT, --min_frac_final FLOAT
                        final minimum fraction of supporting reads per strand (default: 0.76)
  -O INT, --min_hamming_dist INT
                        minimum hamming distance between k-mer starting at the fragment position and any other k-mer in the genome (default: 0, i.e. nothing filtered)
  -P FLOAT, --max_cov_percentile FLOAT
                        maximum sequencing coverage percentile of the variant position (default: 1, i.e. nothing filtered)
  -Q, --require_overlap
                        whether to require support from both forward and reverse mapping reads, i.e. variant within the read1-read2 overlap (default: do not require)
```

### duplex_strand_swapper.py
```
usage: duplex_strand_swapper.py [-h] -i FILES [FILES ...] -o FILE [-u TAG] [-t PREFIX] [-@ INT]

Makes a randomly strand-swapped NanoSeq library from several input NanoSeq libraries. For each fragment start and end position in the input libraries, one set of top strand PCR duplicates (F1R2) and one set of bottom strand PCR
duplicates (F2R1) is randomly selected from two different libraries. If this isn't possible, nothing is output for that fragment. The UMI tag will be overwritten as {index of top strand BAM}:{top strand UMI}_{index of bot strand
BAM}:{bot strand UMI}.

The point of this script is to make all the true mutations uncallable, while keeping the false positive rate the same. Because FPs appear in strands independently, it shouldn't matter how we rearrange the strands, the likelihood of
pairing two strands with the same FP remains the same.

Example:
For one fragment start and end position, three libraries have 2 fragments with reads from both the top and bottom strand. Out of the 6 fragments, we randomly select sample_2 umi_2 to contribute the top strand reads. Out of the 4
fragments not in sample_2, we randomly select sample_0 umi_2 to contribute the bottom strand reads. One fragment will be output with the reads of sample_2 umi_2's top strand and sample_0 umi_2's bottom strand, and the umi will be
"2:umi_2_0:umi_2".

options:
  -h, --help            show this help message and exit
  -i FILES [FILES ...], --input FILES [FILES ...]
                        input NanoSeq BAM files. Must be coordinate sorted and indexed
  -o FILE, --output FILE
                        output strand-swapped BAM file
  -u TAG, --umi TAG     name of tag which contains the fragment umi
  -t PREFIX, --tmp PREFIX
                        prefix to add to temporary files. Will make any directories which do not exist (default: tmp/)
  -@ INT, --threads INT
                        number of threads. Each chromosome can be processed by a separate thread. (default: 1)
```

### duplex_coverage.py
```
usage: duplex_coverage.py [-h] [-t INT] [-s INT] [-q FLOAT] [-@ INT] [-b PREFIX] [-c FLOAT,INT,INT,INT] [-e INT] [-l INT] [-r] [-u UMI] INPUT_FILE OUTPUT_PREFIX

Calculates the callable fragment coverage for duplex-seq. Takes a BAM file of duplex-seq reads as input and determines the number of callable fragments overlapping each position in the genome, in other words, the number of
opportunities a mutation could have been detected and pass all duplex-seq filters at each position. This script is not 100% accurate, as it does not consider indels in reads, and indels alter the reference bases a read spans.

positional arguments:
  INPUT_FILE            BAM file to call from. Must be coordinate sorted and have an index file
  OUTPUT_PREFIX         prefix to append to output numpy arrays, final outputs will be {prefix}{chrom}.npy for each chromosome

options:
  -h, --help            show this help message and exit
  -t INT, --min_total_cov INT
                        number of covering reads required to call variants within a fragment (default: 6)
  -s INT, --min_orig_cov INT
                        number of reads required from each original strand to call variants within a fragment (default: 2)
  -q FLOAT, --min_mq FLOAT
                        minimum average mq of pcr duplicate reads required to call variants within a fragment (default: 10)
  -@ INT, --threads INT
                        chromosomes are divided between threads. Cannot use more threads than num chromosomes + 1 (default: 1)
  -b PREFIX, --blacklist PREFIX
                        prefix of the blacklist files generated by duplex_blacklist_regions. Files should be "{prefix}{contig}_{type}.npy" where type is coverage_percentile, duplicate_distance, poly_at_length, and poly_repeat_length
                        (default: do not use blacklist)
  -c FLOAT,INT,INT,INT, --blacklist_cutoffs FLOAT,INT,INT,INT
                        cutoffs to use for the blacklist. Same order as listed in "-b" help section. Values equal to the cutoff are counted as callable (default: 0.007,1,-1,7,4)
  -e INT, --end_dist INT
                        minimum distance from the fragment end required to call a mutation. The end mismatch rate filter makes fragment ends uncallable (default: 5)
  -l INT, --max_frag_len INT
                        maximum fragment length before a fragment is uncallable (default: 300)
  -r, --require_overlap
                        whether read1-read2 overlap is required to call a variant. WARNING: setting this only gives the correct output when reads are 150bp (default: not required)
  -u UMI, --umi UMI     Consider the UMI in the provided BAM tag when grouping pre-PCR fragments (default: no UMI)
```

### duplex_metadata.py
```
usage: duplex_metadata.py [-h] -i BAM [BAM ...] -o PREFIX [-u TAG] [-r SEPARATOR,INDEX] [-g CHROM:START-END] [-s FLOAT [FLOAT ...]] [-B INT] [-C INT] [-N FLOAT] [-@ INT] [--tmp PREFIX]

Outputs metadata statistics and figures on a list of NanoSeq BAM files. Outputs {prefix}metadata.tsv with columns:
sample rep sample_read_pairs_whole_genome sample_conflict_rate rep_read_pairs rep_frags rep_frags_both_strands rep_frags_call_overlap rep_frags_call_whole - sample: name of BAM file
- rep: name of replicate extracted from the UMI tag (optional)
- sample_read_pairs_whole_genome: total number of concordant read pairs in the sample. All other statistics are for only the region specified with --region
- sample_conflict_rate: the estimated fraction of fragments in the sample which will not be callable due to another fragment at the same alignment position with sufficient reads
- rep_read_pairs: concordant read pairs in the replicate (in --region)
- rep_frags: number of fragments in the replicate (in --region)
- rep_frags_both_strands: number of fragments in the replicate (in --region) with both strands sequenced
- rep_frags_call_overlap: number of fragments in the replicate (in --region) with enough coverage to call in the read1-read2 overlap (>=B/2 read pairs per strand and >=C/2 read pairs total)
- rep_frags_call_whole: number of fragments in the replicate (in --region) with enough coverage to call the whole fragment (>=B read pairs per strand and >=C read pairs total)

Outputs the following figures:
- {prefix}conflict_rates.svg: scatterplot showing the relationship between fragment count and conflict rate as well as the line used to impute conflict rate for samples with one replicate
- {prefix}subsampling_efficiency.svg: used to estimate dilution efficiency. Replicates libraries have their reads subsampled to various levels and sequencing efficiency (callable fragments / aligned read pair) is calculated for each
level. An optimally diluted library will plateau at 1, while an oversequenced library will peak at a subsample level < 1. Note: callable fragments here does not take into account mapping quality or blacklists.
- {prefix}reads_per_molecule.svg: histogram of the number of reads per molecule.

Lastly, a dictionary of the conflict rates is written to stdout for copying into the duplex_mutation_figures notebook.

options:
  -h, --help            show this help message and exit
  -i BAM [BAM ...], --input BAM [BAM ...]
                        input BAM files of duplex-seq libraries
  -o PREFIX, --output PREFIX
                        prefix to prepend the output metadata table and figures
  -u TAG, --umi TAG     BAM tag containing the UMI (default: no UMI)
  -r SEPARATOR,INDEX, --replicate SEPARATOR,INDEX
                        should be set if a single BAM file contains multiple replicate samples. These replicates will be treated separately in the metadata output and are necessary for calculating an estimated conflict rate for a BAM
                        file. The replicate name is extracted from the UMI tag by splitting on SEPARATOR and taking the substring at INDEX. e.g. with "--umi RX --sample_extract _,1", an RX tag of "CACGTC_library1" will extract
                        "library1" as the replicate name. Pass an empty SEPARATOR to indicate that the UMI is the replicate name (default: treat each BAM as one replicate)
  -g CHROM:START-END, --region CHROM:START-END
                        consider only reads with start positions within this region. Significantly reduces runtime (default: use whole genome)
  -s FLOAT [FLOAT ...], --subsample_levels FLOAT [FLOAT ...]
                        subsample the library to these fractions of reads when visualizing dilution efficiency (default: 0.2 0.4 0.6 0.8 1)
  -B INT, --min_strand_cov INT
                        filter cutoff for calling variants, used for determining fragments with enough coverage to call: minimum reads per strand overlapping variant (default: 2)
  -C INT, --min_total_cov INT
                        filter cutoff for calling variants, used for determining fragments with enough coverage to call: minimum total reads overalpping variant (default: 6)
  -N FLOAT, --min_frac_final FLOAT
                        filter cutoff for calling variants, used in calculating conflict rate: final minimum fraction of supporting reads per strand (default: 0.76)
  -@ INT, --threads INT
                        multiprocessing threads for parallel processing of BAM files (default: 1)
  --tmp PREFIX          prefix for temporary files (default: ./)
```

### visualize_duplex_filters.py
```
usage: visualize_duplex_filters.py [-h] -v TSV1 [TSV1 ...] -V TSV1 [TSV1 ...] -w PREFIX1 [PREFIX1 ...] -W PREFIX1 [PREFIX1 ...] -o FILE [-A FLOAT,FLOAT:FLOAT:FLOAT] [-B INT,INT:INT:INT] [-C INT,INT:INT:INT] [-D INT,INT:INT:INT]
                                   [-E FLOAT,FLOAT:FLOAT:FLOAT] [-F FLOAT,FLOAT:FLOAT:FLOAT] [-G INT,INT:INT:INT] [-H FLOAT,FLOAT:FLOAT:FLOAT] [-I INT,INT:INT:INT] [-J INT,INT:INT:INT] [-K INT,INT:INT:INT] [-L INT,INT:INT:INT]
                                   [-M INT,INT:INT:INT] [-N FLOAT,FLOAT:FLOAT:FLOAT] [-O INT,INT:INT:INT] [-P FLOAT,FLOAT:FLOAT:FLOAT] [-Q] [-k]

Visualizes how adjusting the filter thresholds of filter_duplex_variants.py will affect precision and recall. Uses rate of C>T mutations and mutations called in "swapped" control samples to estimate the false positive rate and number
of true positives passing filters. Default values and tested range of filters can be modified by passing values DEFAULT,START:END:STEP to the filter arguments. END is non-inclusive. False positive rate of a given set of filters is
calculated as (passing variants in swapped samples / callable coverage of swapped samples) / (passing variants in real samples / callable coverage of real samples). Note: the callable coverage is static and does not adapt to changes
in filters, so it's best to use swapped samples made from the real samples to minimize any error introduced by this. True positive count is estimated as the number of passing variants in real samples * (1 - false positive rate). C>T
rate refers to the fraction of passing variants in real samples which are C>T or G>A.|n Note: this script may require a large amount of memory, as it needs to hold all the input TSVs in memory. A preliminary filter K is applied to
reduce memory usage, so setting the END value of -K lower will reduce memory footprint.

options:
  -h, --help            show this help message and exit
  -v TSV1 [TSV1 ...], --real_vars TSV1 [TSV1 ...]
                        variant TSV(s) of "real" (non-swapped) samples output by add_duplex_filter_columns.py
  -V TSV1 [TSV1 ...], --swapped_vars TSV1 [TSV1 ...]
                        variant TSV(s) output running add_duplex_filter_columns.py on "swapped" samples generated by duplex_strand_swapper.py
  -w PREFIX1 [PREFIX1 ...], --real_cov PREFIX1 [PREFIX1 ...]
                        prefix of coverage arrays for the same files as -v, generated by duplex_coverage.py
  -W PREFIX1 [PREFIX1 ...], --swapped_cov PREFIX1 [PREFIX1 ...]
                        prefix of coverage arrays for the same files as -V, generated by duplex_coverage.py
  -o FILE, --output FILE
                        name of the output figure (file format is inferred)
  -A FLOAT,FLOAT:FLOAT:FLOAT, --min_frac_init FLOAT,FLOAT:FLOAT:FLOAT
                        initial minimum fraction of supporting reads per strand (default: 0.24,0:1:0.05)
  -B INT,INT:INT:INT, --min_strand_cov INT,INT:INT:INT
                        minimum reads per strand overlapping variant (default: 2,0:10:1)
  -C INT,INT:INT:INT, --min_total_cov INT,INT:INT:INT
                        minimum total reads overalpping variant (default: 6,0:20:1)
  -D INT,INT:INT:INT, --min_high_bq INT,INT:INT:INT
                        minimum supporting reads with high BQ (>=30) (default: 4,0:20:1)
  -E FLOAT,FLOAT:FLOAT:FLOAT, --min_avg_mq FLOAT,FLOAT:FLOAT:FLOAT
                        minimum average MQ of supporting reads (default: 10,0:44:1)
  -F FLOAT,FLOAT:FLOAT:FLOAT, --max_end_mismatch FLOAT,FLOAT:FLOAT:FLOAT
                        maximum mismatch rate between a position and a fragment end before clipping (default: 0.21,0:1:0.02)
  -G INT,INT:INT:INT, --min_indel_end_dist INT,INT:INT:INT
                        minimum distance between indels and a fragment end (default: 6,0:20:1)
  -H FLOAT,FLOAT:FLOAT:FLOAT, --min_cov_percentile FLOAT,FLOAT:FLOAT:FLOAT
                        minimum sequencing coverage percentile of the variant position (default: 0.007,0:0.2:0.01)
  -I INT,INT:INT:INT, --max_poly_at INT,INT:INT:INT
                        maximum length of a poly-A/T repeat ending at the variant position (default: 7,0:10:1)
  -J INT,INT:INT:INT, --max_poly_rep INT,INT:INT:INT
                        maximum length of a poly-di/trinucleotide repeat end at the variant position (default: 4,0:10:1)
  -K INT,INT:INT:INT, --max_control_sup INT,INT:INT:INT
                        maximum number of control sample fragments supporting the variant. The END value of this argument affects how much memory the script requires with lower values reducing memory footprint (default: 3,0:40:1)
  -L INT,INT:INT:INT, --max_frag_len INT,INT:INT:INT
                        maximum fragment length (default: 300,100:2000:100)
  -M INT,INT:INT:INT, --max_contamination_dist INT,INT:INT:INT
                        maximum distance between two variants passing prior filters (A-L) within a single fragment (default: 4,0:10:1)
  -N FLOAT,FLOAT:FLOAT:FLOAT, --min_frac_final FLOAT,FLOAT:FLOAT:FLOAT
                        final minimum fraction of supporting reads per strand (default: 0.76,0:1:0.05)
  -O INT,INT:INT:INT, --min_hamming_dist INT,INT:INT:INT
                        minimum hamming distance between k-mer starting at the fragment position and any other k-mer in the genome (default: 0,0:10:1)
  -P FLOAT,FLOAT:FLOAT:FLOAT, --max_cov_percentile FLOAT,FLOAT:FLOAT:FLOAT
                        maximum sequencing coverage percentile of the variant position (default: 1, i.e. nothing filtered)
  -Q, --require_overlap
                        whether to require support from both forward and reverse mapping reads, i.e. variant within the read1-read2 overlap. Unlike the other filters, this one is just a flag which determines whether to require overlap
                        while varying the other filters (default: do not require)
  -k, --keep_dups       by default, filtered mutations observed in mutliple fragments within a TSV (-v or -V) will only count as a single mutation. Setting this flag will instead count each observation of the mutation. Setting this
                        flag may be useful if TSVs were merged before being passed to -v or -V, and duplicates could represent real mutations or false positives. Note: if the same mutation is found in three different TSVs, it will
                        always count as three mutations.
```

### duplex_subsampler.py
```
usage: duplex_subsampler.py [-h] -o FILE [-s INT] [-t INT] [-c] [-r FLOAT] [-u TAG] [--tmp PREFIX] [-@ INT] FILE

Subsamples fragments/reads from a NanoSeq library BAM.

positional arguments:
  FILE                  input NanoSeq BAM file. Must be coordinate sorted and indexed

options:
  -h, --help            show this help message and exit
  -o FILE, --output FILE
                        subsampled BAM file
  -s INT, --required_strand_coverage INT
                        number of reads in each strand required to output the fragment (default: 0)
  -t INT, --required_total_coverage INT
                        number of total reads in the fragment required to output (default: 0)
  -c, --cap             caps the number of reads output for a fragment at the minimum required to pass the requirements reads to be output are selected at random without replacement (default: output all reads for passing fragments)
  -r FLOAT, --random_frac FLOAT
                        fraction of fragments passing requirements to output (default: 1 (output all passing fragments))
  -u TAG, --umi TAG     Consider the UMI in the provided BAM tag when grouping fragments as PCR dups (default: no UMI)
  --tmp PREFIX          prefix to add to temporary files. Will make any directories which do not exist (default: tmp/)
  -@ INT, --threads INT
                        number of threads. Each chromosome can be processed by a separate thread. (default: 1)
```
