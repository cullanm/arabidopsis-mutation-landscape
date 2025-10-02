sys.stderr.write('loading align.smk\n')
import pandas as pd
import glob
import os
import random

shell('''
mkdir -p data/{{align,multiqc,ref,fastqc,metadata,methyl}}
mkdir -p logs
mkdir -p tmp
mkdir -p benchmarks
''')

# note: wildcards constraints are not shared by included files
wildcard_constraints:
	srr = '([SDE]RR[0-9]+)',	# NCBI SRR accessions
	read = '((_1)|(_2))?',		# strand suffixes for pair end reads (Song 2016 dataset is all single end)
	trim = '(_trim)?',			# '_trim' or nothing
	group = '[^_]*',			# no _ allowed
	subgroup = '[^_]*',			# no _ allowed
	rep = '[^_]*'				# no _ allowed

# gunzip a raw fastq file, only used for local data, not SRR downloads
rule gunzip_fastq:
        input: 'data/fastq/{group}_{subgroup}_{rep}{read}.fastq.gz'
        output: 'data/fastq/{group}_{subgroup}_{rep}{read}.fastq'
        log: 'logs/gunzip_fastq_{group}_{subgroup}_{rep}{read}.out'
        shell: 'gunzip {input}' # pigz can't multithread decompress

# download reference fasta from url in config file
rule fetch_ref_fasta:
	output: 'data/ref/ref.fa'
	log: 'logs/fetch_ref_fasta.out'
	params: gunzip = lambda w: 'true' if config['ref_fasta'][-3:] == '.gz' else 'false'
	resources:
		downloads = 1
	shell: # download and fix chromosome names to match GFF
		# raw string ignores escape chars
		r'''
		wget --no-verbose -O {output}.gz {config[ref_fasta]} &> {log}
		if {params.gunzip}; then
			gunzip {output}.gz
		else
			mv {output}.gz {output}
		fi
		if {config[do_tair_chr_rename]}; then
			sed -i 's/>\([0-9]\)/>Chr\1/g' {output}
			sed -i 's/>chloroplast/>ChrC/g' {output}
			sed -i 's/>M/>ChrM/g' {output}
		fi
		'''

# download reference annotation (gff) from url in config file
rule fetch_ref_annotation:
	output: 'data/ref/ref.' + config['ref_annotation_file_ext']
	log: 'logs/fetch_ref_annotation.out'
	params: 
		url = lambda w: config['ref_annotation'],
		gunzip = lambda w: 'true' if config['ref_annotation'][-3:] == '.gz' else 'false'
	resources:
		downloads = 1
	shell: 
		'''
		wget --no-verbose -O {output}.gz {params.url} &> {log}
		if {params.gunzip}; then
			gunzip {output}.gz
		else
			mv {output}.gz {output}
		fi
		'''
	
# generate a bowtie2 index from the reference fasta
rule ref_bowtie_index:
	input: 'data/ref/ref.fa'
	output: directory('data/ref/bowtie_index')
	log: 'logs/ref_bowtie_index.out'
	threads: 16
	shell: # have tried making a faster index by changing index to have a higher memory footprint, but it made alignment slightly slower
		'''
		rm -r {output} &> null || true
		mkdir -p {output} &> null || true
		bowtie2-build --threads {threads} {input} {output}/index &> {log}
		'''

# generate a star index from reference fasta and annotation
rule ref_star_index:
	input: 
		fasta = 'data/ref/ref.fa',
		gtf = 'data/ref/ref.' + config['ref_annotation_file_ext']
	output: directory('data/ref/star_index')
	log: 'logs/ref_star_index.out'
	params:
		ann_correction = '--sjdbGTFtagExonParentTranscript Parent' if config['ref_annotation_file_ext'] in ['gff3', 'gff'] else ''
	threads: 16
	shell:
		'''
		rm -r {output} &> null || true
		mkdir -p {output} &> null || true
		STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} \
			--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} {params.ann_correction} &> {log}
		'''

# generate a methylpy index from reference fasta
rule ref_methylpy_index:
	input: 'data/ref/ref.fa'
	output: 
		f = directory('data/ref/methylpy_index_f'),
		r = directory('data/ref/methylpy_index_r'),
	log: 'logs/ref_methylpy_index.out'
	params: prefix = lambda w, output: output.f[:-2]
	threads: 16
	shell: 
		'''
		rm -r {output.f} {output.r} &> null || true
		mkdir -p {output.f} {output.r} &> null || true
		methylpy build-reference --num-procs {threads} --input-files {input} --output-prefix {params.prefix} &> {log}
		'''

# prefetch an srr using SRA-tools
rule fetch_srr:
	output: temp(directory('data/fastq/{srr}'))
	log: 'logs/fetch_srr_{srr}.out'
	priority: -1	# don't fetch more if there's anything else we can do
	resources:
		downloads = 1	# restrict the number of downloads the job can be performing at once
	shell: 'prefetch {wildcards.srr} --max-size 100GB -O data/fastq &> {log}'

# extract the fastq file from a single end prefetched SRR using SRA-tools
rule extract_srr_single:
	input: 'data/fastq/{srr}/'
	output: temp('data/fastq/{srr}.fastq')
	log: 'logs/extract_srr_single_{srr}.out'
	threads: 4
	shell: 'fasterq-dump {input} -O data/fastq -t tmp/ -e {threads} &> {log}'

# extract the forward and reverse fastq files from a paired end prefetched SRR using SRA-tools
rule extract_srr_paired:
	input: 'data/fastq/{srr}/'
	output:
		f = temp('data/fastq/{srr}_1.fastq'),		# forward
		r = temp('data/fastq/{srr}_2.fastq')		# reverse
	log: 'logs/extract_srr_paired_{srr}.out'
	threads: 4
	shell: 'fasterq-dump {input} -O data/fastq -t tmp/ -e {threads} &> {log}'

# return the name of an SRR fastq file
def find_srr_fastq(w):
	r = find_sample(w)
	return 'data/fastq/{}{}.fastq'.format(r.run, w.read)

# rename an srr fastq from its accession to something more readable
rule rename_srr:
	input: find_srr_fastq
	output: 'data/fastq/{group}_{subgroup}_{rep}{read}.fastq'
	shell: 'mv {input} {output}'

# extract umi from a pair of fastq files
rule extract_umi:
	input: 
		f = 'data/fastq/{group}_{subgroup}_{rep}_1.fastq',
		r = 'data/fastq/{group}_{subgroup}_{rep}_2.fastq'
	output: 
		f = 'data/fastq/{group}_{subgroup}_{rep}_1_extracted.fastq',
		r = 'data/fastq/{group}_{subgroup}_{rep}_2_extracted.fastq'
	log: 'logs/extract_umi_{group}_{subgroup}_{rep}.out'
	shell: 'python {script_dir}/extract_umi.py {config[extract_umi_args]} {input.f} {input.r} {output.f} {output.r} &> {log}'

# determine whether the normal or umi extracted fastqs are needed for a sample
# def find_extracted_fastqs(w):
	# r = find_sample(w)
	# if r.umi == True:
		# return [f'data/fastq/{r.group}_{r.subgroup}_{r.rep}_1_extracted.fastq', f'data/fastq/{r.group}_{r.subgroup}_{r.rep}_2_extracted.fastq']
	# else:
		# return [f'data/fastq/{r.group}_{r.subgroup}_{r.rep}_1.fastq', f'data/fastq/{r.group}_{r.subgroup}_{r.rep}_2.fastq']

# trim fastq with fastp
rule trim_reads_paired:
	input: 
		R1 = lambda w: f'data/fastq/{r.group}_{r.subgroup}_{r.rep}_1{"_extracted" if find_sample(w).umi else ""}.fastq',
		R2 = lambda w: f'data/fastq/{r.group}_{r.subgroup}_{r.rep}_2{"_extracted" if find_sample(w).umi else ""}.fastq'
	output:
		R1 = temp('data/fastq/{group}_{subgroup}_{rep}_1_trim.fastq'),
		R2 = temp('data/fastq/{group}_{subgroup}_{rep}_2_trim.fastq')
	log: 
		log = 'logs/trim_reads_paired_{group}_{subgroup}_{rep}.html',
		out = 'logs/trim_reads_paired_{group}_{subgroup}_{rep}.out'
	threads: 16
	params:
		adapters = lambda w: '--adapter_sequence={} --adapter_sequence_r2={}'.format(config['trim_adapter_1'], config['trim_adapter_2']) if len(config['trim_adapter_1']) > 0 else ''
	benchmark: 'benchmarks/trim_reads_paired_{group}_{subgroup}_{rep}.tsv'
	shell: 'fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} -h {log.log} -w {threads} {params.adapters} &> {log.out}'

# trim fastq with fastp
rule trim_reads_single:
	input: 'data/fastq/{group}_{subgroup}_{rep}.fastq'
	output: temp('data/fastq/{group}_{subgroup}_{rep}_trim.fastq')
	log: 
		log = 'logs/trim_reads_single_{group}_{subgroup}_{rep}.html',
		out = 'logs/trim_reads_single_{group}_{subgroup}_{rep}.out'
	threads: 8
	benchmark: 'benchmarks/trim_reads_single_{group}_{subgroup}_{rep}.tsv'
	shell: 'fastp -i {input} -o {output} -h {log.log} -w {threads} &> {log.out}'

# run fastqc on a fastq file
rule fastqc:
	input: 'data/fastq/{group}_{subgroup}_{rep}{read}{trim}.fastq'
	output: directory('data/fastqc/{group}_{subgroup}_{rep}{read}{trim}')
	log: 'logs/fastqc_{group}_{subgroup}_{rep}{read}{trim}.out'
	threads: 8
	benchmark: 'benchmarks/fastqc_{group}_{subgroup}_{rep}{read}{trim}.tsv'
	shell: 
		'''
		mkdir -p {output}
		fastqc --dir tmp/ -o {output} -t {threads} {input} &> {log}
		'''
		
# return all pre and post trim fastqc output directories
def find_all_fastqcs(w):
	ret = []
	for i, r in df_samples.iterrows():
		if r.paired:
			ret.append('data/fastqc/{}_{}_{}_1/'.format(r.group, r.subgroup, r.rep))
			ret.append('data/fastqc/{}_{}_{}_2/'.format(r.group, r.subgroup, r.rep))
			ret.append('data/fastqc/{}_{}_{}_1_trim/'.format(r.group, r.subgroup, r.rep))
			ret.append('data/fastqc/{}_{}_{}_2_trim/'.format(r.group, r.subgroup, r.rep))
		else:
			ret.append('data/fastqc/{}_{}_{}/'.format(r.group, r.subgroup, r.rep))
			ret.append('data/fastqc/{}_{}_{}_trim/'.format(r.group, r.subgroup, r.rep))
	return ret

# merge all fastqc outputs into one
rule multiqc:
	input: find_all_fastqcs
	output: 'data/multiqc/multiqc_report.html'
	log: 'logs/multiqc.out'
	params: 
		in_dir = lambda w, input: input[0].rsplit('/', 1)[0],	# won't work if fastqcs are not all in the same folder
		out_dir = lambda w, output: output[0].rsplit('/', 1)[0]
	conda: 'envs/multiqc.yaml'
	resources: mem_mb=6144
	shell: 'multiqc -f --interactive -o {params.out_dir} {params.in_dir} &> {log}'		# -f means overwrite previous reports

# return trimmed fastq name(s) given group, subgroup, and rep in wildcards
def find_trimmed_fastq(w):
	r = find_sample(w)
	if r.paired == True:
		return ['data/fastq/{}_{}_{}_1_trim.fastq'.format(w.group, w.subgroup, w.rep), 
				'data/fastq/{}_{}_{}_2_trim.fastq'.format(w.group, w.subgroup, w.rep)]
	else:
		return ['data/fastq/{}_{}_{}_trim.fastq'.format(w.group, w.subgroup, w.rep)]

# align trimmed fastq(s) to genome with bowtie2
rule align_fastq_bowtie:
	input:
		reads = find_trimmed_fastq,
		ref_index =  'data/ref/bowtie_index/'
	output: temp('data/align/{group}_{subgroup}_{rep}_bowtie.sam')
	log: 'logs/align_fastq_bowtie_{group}_{subgroup}_{rep}.out'
	params:
		reads = lambda w, input: '-1 {} -2 {}'.format(*input.reads) if len(input.reads) == 2 else '-U {}'.format(*input.reads)
	threads: 16
	benchmark: 'benchmarks/align_fastq_bowtie_{group}_{subgroup}_{rep}.tsv'
	shell: 'bowtie2 -x {input.ref_index}/index {params.reads} -S {output} -p {threads} {config[bowtie_args]} &> {log}'

# align trimmed fastq(s) to genome with star
rule align_fastq_star:
	input:
		reads = find_trimmed_fastq,
		ref_index =  'data/ref/star_index/'
	output: 'data/align/{group}_{subgroup}_{rep}_star.sam'
	log: 
		out = 'logs/align_fastq_star_{group}_{subgroup}_{rep}.out',
		log = 'logs/align_fastq_star_{group}_{subgroup}_{rep}.log',
		progress = 'logs/align_fastq_star_{group}_{subgroup}_{rep}_progress.log',
		summary = 'logs/align_fastq_star_{group}_{subgroup}_{rep}_final.log',
		table = 'logs/align_fastq_star_{group}_{subgroup}_{rep}_SJ.tab'
	params: prefix = 'data/align/{group}_{subgroup}_{rep}_'
	threads: 16
	benchmark: 'benchmarks/align_fastq_star_{group}_{subgroup}_{rep}.tsv'
	shell: # this works for paired and single, as lists are converted to " ".join()
		'''
		STAR --runThreadN {threads} --genomeDir {input.ref_index} --readFilesIn {input.reads} \
			--outFileNamePrefix {params.prefix} {config[star_args]} &> {log}
		mv {params.prefix}Aligned.out.sam {output}
		mv {params.prefix}Log.out {log.log}
		mv {params.prefix}Log.progress.out {log.progress}
		mv {params.prefix}Log.final.out {log.summary}
		mv {params.prefix}SJ.out.tab {log.table}
		'''

# align trimmed fastq(s) to genome with methylpy
rule align_fastq_methylpy:
	input:
		reads = find_trimmed_fastq,
		ref_f = 'data/ref/methylpy_index_f/',
		ref_r = 'data/ref/methylpy_index_r',
		ref = 'data/ref/ref.fa'
	output: 
		bam = 'data/align/{group}_{subgroup}_{rep}_methylpy.bam',
		allc = 'data/methyl/{group}_{subgroup}_{rep}.allc'
	log: 'logs/align_fastq_methylpy_{group}_{subgroup}_{rep}.out'
	params: 
		out_path = lambda w, output: output.bam.rsplit('/', 1)[0],
		reads = lambda w, input: 'paired-end-pipeline --read1-files {} --read2-files {}'.format(*input.reads) if len(input.reads) == 2 else 'single-end-pipeline --read-files {}'.format(*input.reads)
	threads: 16
	benchmark: 'benchmarks/align_fastq_methylpy_{group}_{subgroup}_{rep}.tsv'
	shell: 
		'''
		methylpy {params.reads} --sample {wildcards.group}_{wildcards.subgroup}_{wildcards.rep} \
		--forward-ref {input.ref_f} --reverse-ref {input.ref_r} --ref-fasta {input.ref} \
		--path-to-output {params.out_path} --num-procs {threads} \
		--trim-reads False --remove-chr-prefix False --compress-output False --keep-temp-files True \
		{config[methylpy_args]} &> {log}
		mv {params.out_path}/{wildcards.group}_{wildcards.subgroup}_{wildcards.rep}_processed_reads.bam {output.bam}
		mv {params.out_path}/allc_{wildcards.group}_{wildcards.subgroup}_{wildcards.rep}.tsv {output.allc}
		'''
# command for converting an allc to a bed: awk -F '\t' '{print $1 "\t" $2 - 1 "\t" $2 "\t" substr($0,index($0,$3))}'
# allc files are 1 based, BEDs are 0 based

# request the unfiltered sam aligned by star or bowtie2 depending on sample_table
def find_raw_sam(w):
	r = find_sample(w)
	suffix = 'bam' if r.aligner == 'methylpy' else 'sam'
	return f'data/align/{r.group}_{r.subgroup}_{r.rep}_{r.aligner}.{suffix}'

# group mates together with samtools collate
# rule sam_collate:
    # input: find_raw_sam
    # output: temp('data/align/{group}_{subgroup}_{rep}_collate.sam')
    # log: 'logs/sam_collate_{group}_{subgroup}_{rep}.out'
    # threads: 8 # this rule ran out of memory on larger samples
    # wildcard_constraints: rep = '(?!merged).*'	# don't allow anything that starts with "merged", as that leads to ambiguous rule. The real regex is '^(?!merged).*$' but snakemake doesn't like anchors
    # shell: 'samtools collate -@ {threads} -o {output} {input} &> {log}'

# add mate info with samtools fixmate (needed for samtools markdup)
rule sam_fixmate:
    # input: 'data/align/{group}_{subgroup}_{rep}_collate.sam'
    input: find_raw_sam # Bowtie2 outputs in collated order, so the collate rule is unneccessary. not sure if STAR and methylpy do the same though
    output: temp('data/align/{group}_{subgroup}_{rep}_fixmate.bam')
    log: 'logs/sam_fixmate_{group}_{subgroup}_{rep}.out'
    threads: 8
    wildcard_constraints: rep = '(?!merged).*'
    shell: 'samtools fixmate -@ {threads} -u -m {input} {output} &> {log}'

# position sort with sambamba
rule sam_sort:
    input: 'data/align/{group}_{subgroup}_{rep}_fixmate.bam'
    output: temp('data/align/{group}_{subgroup}_{rep}_sort.bam')
    log: 'logs/sam_sort_{group}_{subgroup}_{rep}.out'
    threads: 12  # I think the -m option is memory per thread, as setting it over 1G can run out of memory
    resources: mem_mb=24576
    shell: 'sambamba sort -t {threads} -m 512MiB --compression-level 1 -o {output} {input} &> {log}'

# if the sample has a UMI, add the read group tag ("UG") and corrected UMI ("RX") + the replicate name, requires a coordinate sorted sam
rule sam_group_umi:
	input: 'data/align/{group}_{subgroup}_{rep}_sort.bam'
	output: temp('data/align/{group}_{subgroup}_{rep}_umigroup.bam')
	log: 'logs/sam_group_umi_{group}_{subgroup}_{rep}.out'
	# params: umi = lambda w: check_sample_bool(w, 'umi')
	threads: 8
	shell: 'python {script_dir}/group_umis.py -@ {threads} --append _{wildcards.rep} {input} {output} &> {log}'

# returns the sam file to run markdup on depending on if umis need to be grouped first
# def find_sam_markdup_input(w):
	# r = find_sample(w)
	# suffix = 'umigroup' if r.umi == True else 'sort'
	# return f'data/align/{r.group}_{r.subgroup}_{r.rep}_{suffix}.bam'

# mark duplicate reads, taking into account the UMI if umi==true and marking optical duplicates if dedup_opt==true.
# duplicate reads will have the SQ tag if an optical dup and the LB tag if a PCR dup
rule sam_markdup:
    input: lambda w: f'data/align/{r.group}_{r.subgroup}_{r.rep}_{"umigroup" if find_sample(w).umi else "sort"}.bam'
    output: 'data/align/{group}_{subgroup}_{rep}_markdup.bam'
    log: 'logs/sam_markdup_{group}_{subgroup}_{rep}.out'
    threads: 8
    params: req_umi = lambda w: '--barcode-tag RX' if find_sample(w).umi else ''
    shell: 'samtools markdup -@ {threads} {params.req_umi} -d {config[optical_duplicate_dist]} -s {input} {output} &> {log}'

def make_filter_expression(w):
    exp = 'mapping_quality >= 1'    # MQ filter
    r = find_sample(w)
    if r.paired:
        exp += ' and proper_pair'   # concordant align filter
    if r.dedup_opt and r.dedup_pcr:
        exp += " and [dt] == null"  # optical and pcr duplicate filter
    elif r.dedup_opt:
        exp += " and ([dt] == null or [dt] != 'SQ')" # only optical duplicate filter
    elif r.dedup_pcr:
        exp += " and ([dt] == null or [dt] != 'LB')"  # only PCR duplicate filter
    return exp

rule sam_filter:
    input: 'data/align/{group}_{subgroup}_{rep}_markdup.bam'
    output: 'data/align/{group}_{subgroup}_{rep}_filtered.bam'
    log: 'logs/sam_filter_{group}_{subgroup}_{rep}.out'
    threads: 8
    params: filter = make_filter_expression
    shell: 'sambamba view -t {threads} --format bam -l 3 --filter "{params.filter}" -o {output} {input}'

# index any BAM file
rule index_bam:
	input: '{anything}.bam'
	output: '{anything}.bam.bai'
	params: log = lambda w: 'logs/index_bam_{}.out'.format(w.anything.split('/')[-1])
	threads: 8
	shell: 'samtools index -@ {threads} {input} &> {params.log}'

# return BAM files for all replicates of a tf and treatment specified in wildcards
def find_replicate_bams(w):
	reps = df_samples[(df_samples.group == w.group) & (df_samples.subgroup == w.subgroup)].rep
	return ['data/align/{}_{}_{}_filtered.bam'.format(w.group, w.subgroup, rep) for rep in reps]

# merge BAMs with the same group and subgroup
rule merge_replicate_bams:
	input: 
		bams = find_replicate_bams
	output: 'data/align/{group}_{subgroup}_merged.bam'
	log: 'logs/merge_replicate_bams_{group}_{subgroup}.out'
	threads: 8
	shell: 'samtools merge -@ {threads} -r -o {output} {input} &> {log}'

# FIXME this won't work with aligners other than bowtie2
rule alignment_stats:
	input: expand('logs/align_fastq_{group}_{subgroup}_{rep}.out', zip, group=df_samples.group, subgroup=df_samples.subgroup, rep=df_samples.rep)
	output: 'data/metadata/alignment_stats.tsv'
	run: # I think "run" blocks will also rerun any code not in a rule
		data = []
		# for some reason "input" is empty here, so I use the expansion again
		for i in expand('logs/align_fastq_{group}_{subgroup}_{rep}.out', zip, group=df_samples.group, subgroup=df_samples.subgroup, rep=df_samples.rep):
			row = dict()
			with open(i, 'r') as f:
				row['name'] = i.split('_', 2)[-1].rsplit('.', 1)[0]
				for l in f:
					if 'reads; of these' in l:
						row['total'] = l.strip().split()[0]
					elif ') aligned concordantly 0 times' in l:
						row['zero_concordant'] = l.strip().split()[0]
					elif 'concordantly exactly 1 time' in l:
						row['aligned'] = l.strip().split()[0]
					elif 'concordantly >1 times' in l:
						row['multialigned'] = l.strip().split()[0]
					elif 'discordantly 1 time' in l:
						row['discordant'] = l.strip().split()[0]
			data.append(row)
		pd.DataFrame(data).to_csv(output[0], sep='\t', index=False)
