sys.stderr.write('loading duplex.smk\n')
import pandas as pd
import glob
import os
import random

shell('mkdir -p data/{{variant,coverage}}')

# note: wildcards constraints are not shared by included files
wildcard_constraints:
	srr = '([SDE]RR[0-9]+)',	# NCBI SRR accessions
	read = '((_1)|(_2))?',		# strand suffixes for pair end reads (Song 2016 dataset is all single end)
	trim = '(_trim)?',			# '_trim' or nothing
	group = '[^_]*',			# no _ allowed
	subgroup = '[^_]*',			# no _ allowed
	rep = '[^_]*'				# no _ allowed

if sum((df_samples.group == 'duplex') & (df_samples.subgroup == 'swapped')) > 0:
	sys.stderr.write('ERROR: samples cannot have a group of "duplex" and subgroup of "swapped"\n')
	exit()

# make a df of the control samples
df_controls = df_samples[df_samples.control].drop_duplicates('group subgroup'.split())
if len(df_controls) == 0:
	sys.stderr.write('ERROR: no samples labeled as control\n')

# make a coverage bedgraph from a bam file, used in generate_duplex_blacklists
rule coverage_bg:
	input: 'data/align/{anything}.bam'
	output: temp('data/coverage/{anything}.bg')
	log: 'logs/coverage_bg_{anything}.out'
	shell: 'bedtools genomecov -bg -ibam {input} > {output} 2> {log}'

# make blacklist arrays for add_duplex_filter_columns
rule generate_duplex_blacklists:
	input: 
		bgs = [f'data/coverage/{r.group}_{r.subgroup}_merged.bg' for r in df_controls.itertuples()],
		fasta = 'data/ref/ref.fa'
	output: directory('data/ref/duplex_blacklist/') # note this will not have the trailing "/" when converted to a string
	log: 'logs/generate_duplex_blacklists.out'
	threads: 8
	shell: 'python {script_dir}/generate_duplex_blacklists.py -@ {threads} --fasta {input.fasta} --bedgraphs {input.bgs} --output {output}/ &> {log}'

def get_umi(w):
	if w.group == 'duplex' and w.subgroup == 'swapped':
		return '--umi RX' if df_controls.iloc[0].umi else '--umi RG'
	else:
		return '--umi RX' if find_sample(w).umi else ('--umi RG' if w.rep == 'merged' else '')

# get unfiltered NanoSeq variants
rule duplex_caller:
	input: 
		bam = 'data/align/{group}_{subgroup}_{rep}.bam',
		bai = 'data/align/{group}_{subgroup}_{rep}.bam.bai'
	output: 'data/variant/{group}_{subgroup}_{rep}_duplex.tsv'
	log: 'logs/duplex_caller_{group}_{subgroup}_{rep}.out'
	params:
		umi = get_umi
	threads: 8
	shell: 'python {script_dir}/duplex_caller.py -@ {threads} {params.umi} {input.bam} {output} &> {log}'

# count the number of fragments supporting each variant, only consider fragments with at least 2 reads covering the site. 
# Will be used in filtering the NanoSeq variants and determining their abundance in a sample
rule dup_informed_caller:
	input:
		bam = 'data/align/{group}_{subgroup}_{rep}.bam',
		bai = 'data/align/{group}_{subgroup}_{rep}.bam.bai'
	output: 'data/variant/{group}_{subgroup}_{rep}_informed.tsv'
	log: 'logs/dup_informed_caller_{group}_{subgroup}_{rep}.tsv'
	params: 
		umi = get_umi
	threads: 8
	shell: 'python {script_dir}/dup_informed_caller.py -@ {threads} {params.umi} --duplicate_support 2 {input.bam} {output} &> {log}'

rule add_duplex_filter_columns:
	input: 
		dup = 'data/variant/{group}_{subgroup}_{rep}_duplex.tsv',
		inf = 'data/variant/{group}_{subgroup}_{rep}_informed.tsv',
		controls = lambda w: [f'data/variant/{r.group}_{r.subgroup}_merged_informed.tsv' for r in df_controls.itertuples() if not (r.group == w.group and r.subgroup == w.subgroup)],
		blacklists = 'data/ref/duplex_blacklist/' # this also doesn't add the "/" when converted to text
	output: 'data/variant/{group}_{subgroup}_{rep}_added.tsv'
	log: 'logs/add_duplex_filter_columns_{group}_{subgroup}_{rep}.out'
	shell: 'python {script_dir}/add_duplex_filter_columns.py --input {input.dup} --frequency {input.inf} --controls {input.controls} --blacklists {input.blacklists}/ --output {output} &> {log}'

rule filter_duplex_variants:
	input: 'data/variant/{group}_{subgroup}_{rep}_added.tsv'
	output: 'data/variant/{group}_{subgroup}_{rep}_filtered.tsv'
	log: 
		out = 'logs/filter_duplex_variants_{group}_{subgroup}_{rep}.out',
		err = 'logs/filter_duplex_variants_{group}_{subgroup}_{rep}.err'
	resources: mem_mb=4096
	shell: 'python {script_dir}/filter_duplex_variants.py {config[duplex_filter_args]} --input {input}  --output {output} > {log.out} 2> {log.err}'

rule duplex_coverage:
	input: 
		bam = 'data/align/{group}_{subgroup}_{rep}.bam',
		bai = 'data/align/{group}_{subgroup}_{rep}.bam.bai',
		blacklists = 'data/ref/duplex_blacklist/'
	output: directory('data/coverage/{group}_{subgroup}_{rep}_duplex/')
	log: 'logs/duplex_coverage_{group}_{subgroup}_{rep}.out'
	params: 
		umi = get_umi
	threads: 8
	shell: 'python {script_dir}/duplex_coverage.py -@ {threads} {params.umi} --blacklist {input.blacklists}/ {input.bam} {output}/ &> {log}'

def warn_diff_umi(df):
	if sum(df.umi) != len(df) and sum(df.umi) != 0:
		sys.stderr.write('WARNING: not all samples have the same umi column, this will make duplex_metadata and duplex_filters.svg inaccurate or fail\n')

rule duplex_metadata:
	input: 
		bams = [f'data/align/{r.group}_{r.subgroup}_merged.bam' for r in df_samples.drop_duplicates('group subgroup'.split()).itertuples()],
		bais = [f'data/align/{r.group}_{r.subgroup}_merged.bam.bai' for r in df_samples.drop_duplicates('group subgroup'.split()).itertuples()]
	output: directory('data/metadata/duplex_metadata/')
	log: 
		out = 'logs/duplex_metadata.out',
		err = 'logs/duplex_metadata.err'
	params: 
		warn = warn_diff_umi(df_samples),
		umi = '--umi RX' if df_samples.iloc[0].umi else '--umi RG',
		rep = '--replicate _,1' if df_samples.iloc[0].umi else '--replicate ,0'
	threads: 16
	shell: 'python {script_dir}/duplex_metadata.py -@ {threads} {params.umi} {params.rep} --region Chr1:0-1000000 --input {input.bams} --output {output}/ > {log.out} 2> {log.err}'

rule duplex_strand_swapper:
	input:
		bams = [f'data/align/{r.group}_{r.subgroup}_merged.bam' for r in df_controls.itertuples()],
		bais = [f'data/align/{r.group}_{r.subgroup}_merged.bam.bai' for r in df_controls.itertuples()],
	output: 'data/align/duplex_swapped_{x}.bam'
	log: 'logs/duplex_strand_swapper_{x}.out'
	params: 
		warn = warn_diff_umi(df_controls),
		umi = '--umi RX' if df_controls.iloc[0].umi else '--umi RG'
	threads: 8
	shell: 'python {script_dir}/duplex_strand_swapper.py -@ {threads} -k {config[num_swapped_resamples]} {params.umi} --input {input.bams} --output {output} &> {log}'
	
rule visualize_duplex_filters:
	input:
		real_duplex = [f'data/variant/{r.group}_{r.subgroup}_merged_added.tsv' for r in df_controls.itertuples()],
		real_cov = [f'data/coverage/{r.group}_{r.subgroup}_merged_duplex/' for r in df_controls.itertuples()],
		swapped_duplex = [f'data/variant/duplex_swapped_{i}_added.tsv' for i in range(config['num_swapped_samples'])],
		swapped_cov = [f'data/coverage/duplex_swapped_{i}_duplex/' for i in range(config['num_swapped_samples'])]
	output: 'data/metadata/duplex_filters.svg'
	resources: mem_mb=32768
	params: # need to add slashes to the end of the cov prefixes, as snakemake removes them from input
		real_cov = lambda w, input: '/ '.join(input.real_cov) + '/',
		swapped_cov = lambda w, input: '/ '.join(input.swapped_cov) + '/',
	log: 'logs/visualize_duplex_filters.out'
	shell: 'python {script_dir}/visualize_duplex_filters.py --real_vars {input.real_duplex} --swapped_vars {input.swapped_duplex} --real_cov {params.real_cov} --swapped_cov {params.swapped_cov} --output {output} &> {log}'