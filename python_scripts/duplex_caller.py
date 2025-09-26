#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pysam
import re
from collections import defaultdict
from tqdm import tqdm
import functools
import numpy as np
import argparse
import sys
import multiprocessing
import os.path
import os
import time
import random
import gtools


# In[ ]:


parser = argparse.ArgumentParser(description='''
Calls variants from a coordinate sorted BAM file of duplex sequencing reads. Ignores reads that didn't map as a concordant pair. Output will treat back-to-back SNVs and 
indels longer than 1bp as multiple 1bp variants, these will need to be grouped later. For insertions, this means a CG insertion will have two entries in the output,
a *->C and *->CG, with the second entry listing the support for the inserted G only. Note, the vast majority of variants output by this script will not be true variants,
further filtering steps will be necessary. Output is a tsv with the following columns:|n
    - chrom: chromosome/contig|n
    - pos: 0-based position of the variant. For insertions, this is the site in the reference which the insertion occurs before|n
    - ref: reference sequence, "*" for insertion|n
    - alt: alt sequence, "*" for deletion|n
    - frag_start: start position of the original fragment which contains the variant (two variants at the same position may be on different fragments)|n
    - frag_length: length of the original fragment which contains the variant|n
    - frag_umi: UMI of the original fragment (will be empty if umi argument isn't specified)|n
    - f_read1_sup: number of reads generated from the forward strand of the original fragment which support the variant|n
    - f_read1_cov: number of reads generated from the forward strand of the original fragment which overlap the variant|n
    - r_read1_sup: number of reads generated from the reverse strand of the original fragment which support the variant|n
    - r_read1_cov: number of reads generated from the reverse strand of the original fragment which overlap the variant|n
    - f_sup: number of reads mapping to the forward strand which support the variant|n
    - f_cov: number of reads mapping to the forward strand which overlap the variant|n
    - r_sup: number of reads mapping to the reverse strand which support the variant|n
    - r_cov: number of reads mapping to the reverse strand which overlap the variant|n
    - f_read1_conc: number of read pairs from the forward strand of the original fragment where both read1 and read2 support the variant|n
    - r_read1_conc: number of read pairs from the reverse strand of the original fragment where both read1 and read2 support the variant|n
    - mq: mapping quality of the reads supporting the variant. ascii converted (qual + 33)|n
    - base quality: base call quality of the variant bases. ascii converted (qual + 33)|n
    - end_mismatch_rate: maximum mismatch rate of any window starting from a fragment end and containing the variant. e.g. Removing all variants with 
    end_mismatch_rate > 0.5 would mean ignoring the ends of each fragment until there are more matches than mismatches from the position to the fragment end.
    ''', formatter_class=gtools.EscapeFormatter)
parser.add_argument(dest='input', metavar='INPUT_FILE', type=str, 
                   help='BAM file to call from. Must be coordinate sorted and have an index file. Reads must be paired end with no unpaired reads present.')
parser.add_argument(dest='output', metavar='OUTPUT_FILE', type=str,
                   help='output variant file')
parser.add_argument('-@', '--threads', dest='threads', metavar='THREADS_INT', type=int, default=1,
                   help='chromosomes are divided between threads. Cannot use more threads than num chromosomes + 1 (default: 1)')
parser.add_argument('-t', '--tmp', dest='tmp', metavar='TMP_STRING', type=str, default='./',
                   help='prefix to append temporary output files (default: ./)')
parser.add_argument('-u', '--umi', dest='umi', metavar='TAG', type=str, default=None,
                   help='Consider the UMI in the provided BAM tag when grouping pre-PCR fragments (default: no UMI)')
parser.add_argument('-a', '--all', dest='all', action='store_true',
                   help='when set, outputs all variants rather than only those with support in both original strands')
parser.add_argument('-C', '--chunk_size', dest='chunk_size', metavar='INT', type=int, default=10000,
                   help='''No effect on output. Number of reads to process before flushing buffered information. Increasing
                   may speed up runtime, especially for high-coverage samples. Decreasing lowers memory requirement (default: 10000)''')

try: # run this if in a jupyter notebook
    get_ipython()
    # args = parser.parse_args('tests/test_duplex_caller/duplex_caller_test_case.bam tests/test_duplex_caller/out.tsv -@ 1 --all -t tmp/'.split()) # used for testing
    if os.getcwd()[:8] != '/scratch': # switch to the scratch directory where all the data files are
        os.chdir(f'/scratch/cam02551/{os.getcwd().split("/")[-2]}')
    args = parser.parse_args('-@ 6 --umi RG --all data/align/big_parp2-1_merged.bam tmp/test_duplex_caller.tsv'.split()) # used for testing
except: # run this if in a terminal
    args = parser.parse_args()

if not os.path.isfile(args.input):
    sys.stderr.write('ERROR: input file {} does not exist\n'.format(args.input))
    exit()

aln = pysam.AlignmentFile(args.input, 'rb') # load BAM for doing checks

# assert index can be opened
try:
    aln.check_index()
except:
    sys.stderr.write('ERROR: failed to open index file {}\n'.format(args.input + '.bai'))
    exit()

if args.output and '/' in args.output:
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

if args.tmp and '/' in args.tmp:
    os.makedirs(os.path.dirname(args.tmp), exist_ok=True)

# set threads to num chromosomes + 1 if more were allocated
if args.threads > len(aln.get_index_statistics()) + 1:
    args.threads = len(aln.get_index_statistics()) + 1
    sys.stderr.write('WARNING: specified more threads than number of chromosomes + 1. Reducing threads to {}\n'.format(len(aln.get_index_statistics()) + 1))

sys.stderr.write('Running duplex_caller.py with arguments:\n' + '\n'.join([f'{key}={val}' for key, val in vars(args).items()]) + '\n')


# In[ ]:


# confirm that file is coordinate sorted and contains the MD tag and specified UMI
last = ('', 0)
for i, read in enumerate(aln.fetch()):
    # check if UMI tag is present
    if args.umi is not None:
        try:
            read.get_tag(args.umi)
        except KeyError:
            sys.stderr.write(f'ERROR: BAM file does not contain the specified UMI tag "{args.umi}"\n')
            exit()

    # check if MD tag is present
    if not read.is_unmapped and (read.flag & 2): # if properly mapped
        try:
            read.get_tag('MD')
        except KeyError:
            sys.stderr.write('ERROR: BAM file does not contain the "MD" tag\n')
            exit()

    # check if current read aligned before previous read (out of order)
    nex = (read.reference_name, read.reference_start)
    if nex < last:
        sys.stderr.write('ERROR: BAM file is not coordinate sorted\n')
        exit()
    last = nex
    if i > 1000:
        break


# ### Notes:
# - code requires coordinate sorted BAM
# - code requires the BAM to contain the MD tag to function (see https://samtools.github.io/hts-specs/SAMtags.pdf). Tag follows the regex `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`, where a number indicates a match with the reference, a letter indicates the reference base at a mismatch, and a deletion is specified as a carat followed by one or more letters indicating the aligned reference sequence. This is how Bowtie2 produces the MD field. Insertions are not represented in the MD tag.
# - all position variables and output should be zero-based

# In[ ]:


draw_pbars = True # whether to draw tqdm progress bars


# In[ ]:


# these will be the column names ofthe final output, the temporary outputs of call_variants must match this order
output_columns = 'chrom pos ref alt frag_start frag_len frag_umi f_read1_sup f_read1_cov r_read1_sup r_read1_cov f_sup f_cov r_sup r_cov f_read1_conc r_read1_conc mq bq end_mismatch_rate'.split()


# In[ ]:


# calls variants for a single contig
def call_variants(input_filepath, contig, total_reads, contig_len, output_filepath, process_num, lock):
    aln = pysam.AlignmentFile(input_filepath, 'rb')
    args.chunk_size = 10000

    # key is (pos (0 based), ref, alt, frag_start, frag_length, umi)
    variants = defaultdict(lambda: {'f_read1_sup': 0, # supporting reads in original + strand
                                    'r_read1_sup': 0, # supporting reads in original - strand
                                    'f_sup': 0, # supporting reads on + strand
                                    'r_sup': 0, # supporting reads on - strand
                                    'f_read1_conc': 0, # number of read pairs on original + strand where both support
                                    'r_read1_conc': 0, # number of read pairs on original - strand where both support
                                    'mq': '', # mapping qualities of supporting reads
                                    'bq': '', # base qualities of supporting reads
                                    'f_read1_names': set(), # names of supporting reads on original + strand
                                    'r_read1_names': set(), # names of supporting reads on original - strand
                                    'start_mismatch_rates': [], # maximal mismatch rates from variant toward the fragment start
                                    'end_mismatch_rates': [], # maximal mismatch rates from variant toward the fragment end
                                    'start_coverage': 0, # total bases from variant to read start
                                    'end_coverage': 0}) # total bases from variant to read end
    frag_coverage = dict() # key is (frag_start, frag_length, umi, category) where category is f_read1, r_read1, f, or r; value is array of read coverage over the fragment
    # chrom_coverage = np.zeros(contig_len, dtype=np.short) # read coverage of chromosome (all fragments)
    num_read = 0
    cur_pos = -1
    ready_to_flush = False
    output_file = open(output_filepath, 'w') # temporary output file

    # initialize progress bar
    if draw_pbars:
        with lock:
            pbar = tqdm(total=total_reads, desc=contig, position=process_num, leave=True, mininterval=10)

    for read in aln.fetch(contig=contig):
        num_read += 1

        # skip pairs that don't align as expected
        if read.is_unmapped or read.mate_is_unmapped or not (read.flag & 2):
            continue

        # determine the strand of the fragment based on the alignment of read1
        if read.is_read1:
            frag_strand = 'r_read1' if read.is_reverse else 'f_read1'
        else:
            frag_strand = 'f_read1' if read.is_reverse else 'r_read1'

        # determine the fragment start and length
        frag_start = min(read.reference_start, read.next_reference_start)
        frag_len = abs(read.template_length)

        # split md tag into matches (numbers), mismatches (letters), and deletions (carat + letter(s)). non-matches are always separated by a number
        mds = re.findall('[0-9]+|[A-Z]|\^[A-Z]+', read.get_tag('MD'))

        # change reference N mismatches to matches (e.g. ['3', 'N', '3'] -> ['7'])
        new_mds = []
        was_n = False
        for md in mds:
            if md[0] == 'N':
                new_mds[-1] = str(int(new_mds[-1]) + 1)
                was_n = True
            elif was_n:
                new_mds[-1]  = str(int(new_mds[-1]) + int(md))
                was_n = False
            else:
                new_mds.append(md)
        mds = new_mds

        # change deletions of lenth >1 into multiple 1bp deletions (e.g. ['^AGC'] -> ['^A', '^G', '^C'])
        # new_mds = []
        # for md in mds:
        #     if md[0] == '^':
        #         for i in range(1, len(md)):
        #             new_mds.append('^' + md[i])
        #     else:
        #         new_mds.append(md)
        # mds = new_mds

        cur_cig = 0 # keep track of the sum of cigar element lengths
        # find the insertions in the cigar, add them to mds as "*x", where x is the length of the insertion
        for cig in read.cigar:
            if cig[0] == 1: # if insertion
                cur_md = 0 # keep track of sum of md element lengths
                for i, md in enumerate(mds): # find md entry which the insertion falls within
                    if md.isnumeric(): # if match
                        cur_md += int(md)
                        if cur_cig <= cur_md: # if the insertion falls within this md element
                            first = int(md) - (cur_md - cur_cig)
                            second = int(md) - first
                            insert = []
                            if first != '0':
                                insert.append(str(first))
                            insert.append('*' + str(cig[1]))
                            if second != '0':
                                insert.append(str(second))
                            mds = mds[:i] + insert + mds[i+1:]
                            break
                    elif md[0] == '*': # if insertion
                        cur_md += int(md[1:])
                    elif md[0] == '^': # if deletion
                        cur_md += len(md) - 1
                    else: # if mismatch
                        cur_md += 1

            elif cig[0] != 0 and cig[0] != 2: # if not aligned nor deletion
                with lock:
                    sys.stderr.write('WARNING: treating unexpected cigar char "{}" as match in read:\n{}\n'.format(cig[0], read))
            cur_cig += cig[1]

        # remove spacer zeros (e.g. ['A', '0', 'T'] -> ['A', 'T'])
        mds = [md for md in mds if md != '0']

        # example of how mds should now look: ['A', '3', 'C', 'G', '^TA', '9', '*2', '2']

        total_mismatches = 0
        for md in mds:
            if md.isnumeric():
                continue
            elif md[0] == '*':
                total_mismatches += int(md[1:])
            elif md[0] == '^':
                total_mismatches += len(md) - 1
            else:
                total_mismatches += 1

        # get the UMI tag
        if args.umi is None:
            umi = ''
        else:
            umi = read.get_tag(args.umi) 

        # get the worst variant to fragment edge mismatch rate for each md
        edge_mm_rate = np.zeros(len(mds))
        dist = 0
        past_mismatches = 0
        if read.is_reverse: # for reverse reads, only consider the position to the fragment end
            for i, md in enumerate(mds):
                if md.isnumeric():
                    dist += int(md)
                    continue

                r = (total_mismatches - past_mismatches) / (read.query_length - dist)
                if edge_mm_rate[i] < r:
                        edge_mm_rate[i:] = r

                if md[0] == '^':
                    past_mismatches += len(md) - 1
                elif md[0] == '*':
                    past_mismatches += int(md[1:])
                    dist += int(md[1:])
                else:
                    past_mismatches += 1
                    dist += 1
        else: # for forward reads, only consider the position to the fragment start
            for i in range(len(mds) - 1, -1, -1):
                md = mds[i]
                if md.isnumeric():
                    dist += int(md)
                    continue

                r = (total_mismatches - past_mismatches) / (read.query_length - dist)
                if edge_mm_rate[i] < r:
                        edge_mm_rate[:i + 1] = r

                if md[0] == '^':
                    past_mismatches += len(md) - 1
                elif md[0] == '*':
                    past_mismatches += int(md[1:])
                    dist += int(md[1:])
                else:
                    past_mismatches += 1
                    dist += 1


        # find and add all variants using md
        cur_md = 0 # current position in the read
        total_insertion_len = 0
        for i, md in enumerate(mds):
            if md.isnumeric(): # skip over matches
                cur_md += int(md)
            else:
                if md[0] == '^': # if deletion
                    key_pos = read.reference_start + cur_md - total_insertion_len
                    key_ref = md[1:]
                    key_alt = '*'
                    qualities = '' # no base call qualities to add
                    # don't increment cur_md, as deletions don't count towards the read length
                    total_insertion_len -= len(md) - 1
                elif md[0] == '*': # if insertion
                    key_pos = read.reference_start + cur_md - total_insertion_len
                    key_ref = '*'
                    key_alt = read.query_sequence[cur_md:cur_md + int(md[1:])]
                    qualities = ''
                    for qual in read.query_qualities[cur_md:cur_md+int(md[1:])]:
                        qualities += chr(qual + 33)
                    cur_md += int(md[1:])
                    total_insertion_len += int(md[1:])
                else: # if mismatch
                    key_pos = read.reference_start + cur_md - total_insertion_len
                    key_ref = md
                    key_alt = read.query_sequence[cur_md]
                    qualities = chr(read.query_qualities[cur_md] + 33)
                    cur_md += 1

                key = (key_pos, key_ref, key_alt,
                       frag_start, frag_len, 
                       umi)

                # add variants to dictionary
                if frag_strand == 'f_read1':
                    if read.query_name in variants[key]['f_read1_names']: # if we already found this variant in the mate
                        variants[key]['f_read1_conc'] += 1
                    else:
                        variants[key]['f_read1_names'].add(read.query_name)
                    variants[key]['f_read1_sup'] += 1
                else:
                    if read.query_name in variants[key]['r_read1_names']:
                        variants[key]['r_read1_conc'] += 1
                    else:
                        variants[key]['r_read1_names'].add(read.query_name)
                    variants[key]['r_read1_sup'] += 1
                if read.is_reverse:
                    variants[key]['r_sup'] += 1
                    variants[key]['end_mismatch_rates'].append(edge_mm_rate[i])
                else:
                    variants[key]['f_sup'] += 1
                    variants[key]['start_mismatch_rates'].append(edge_mm_rate[i])
                variants[key]['mq'] += chr(read.mapping_quality + 33)
                variants[key]['bq'] += qualities


        # add coverage of fragment to frag_coverage
        frag = (frag_start, frag_len, umi)
        if frag not in frag_coverage: # if not in dictionary yet, add it
            frag_coverage[frag] = {'f_read1':np.zeros(frag[1], dtype=np.short), 
                                   'r_read1':np.zeros(frag[1], dtype=np.short), 
                                   'f':np.zeros(frag[1], dtype=np.short), 
                                   'r':np.zeros(frag[1], dtype=np.short)}
        if read.is_reverse: # FIXME there may be errors here if forward read goes past the end of the reverse read (haven't seen that yet in concordant bowtie2 reads)
            frag_coverage[frag]['r'][frag_len - read.reference_length:] += 1
            frag_coverage[frag][frag_strand][frag_len - read.reference_length:] += 1
        else:
            frag_coverage[frag]['f'][:read.reference_length] += 1
            frag_coverage[frag][frag_strand][:read.reference_length] += 1

        # periodically flush variants and coverage info out of the dictionaries
        if num_read % args.chunk_size == 0:
            if draw_pbars:
                with lock:
                    pbar.update(args.chunk_size)
            ready_to_flush = True

        if (ready_to_flush and read.reference_start > cur_pos) or num_read == total_reads:
            cur_pos = read.reference_start

            # flush variants
            to_flush = []
            for var in variants: 
                if var[0] >= cur_pos and num_read < total_reads: # don't yet flush variants at or past the current read start
                    continue
                to_flush.append(var)
            to_flush.sort() # this should ensure the output is sorted by position-ref-alt
            for var in to_flush:
                frag = (var[3], var[4], var[5]) # (frag start, frag end, umi)
                pos_in_frag = var[0] - var[3] # variant position - frag start
                f_read1_cov = frag_coverage[frag]['f_read1'][pos_in_frag]
                r_read1_cov = frag_coverage[frag]['r_read1'][pos_in_frag]
                f_cov = frag_coverage[frag]['f'][pos_in_frag]
                r_cov = frag_coverage[frag]['r'][pos_in_frag]
                s_rates = variants[var]['start_mismatch_rates']
                e_rates = variants[var]['end_mismatch_rates']
                avg_start_mm_rate = sum(s_rates) / len(s_rates) if len(s_rates) > 0 else 0
                avg_end_mm_rate = sum(e_rates) / len(e_rates) if len(e_rates) > 0 else 0
                max_mr = '%.3f' % max(avg_start_mm_rate, avg_end_mm_rate)

                v = variants[var]
                if args.all or (v['f_read1_sup'] > 0 and v['r_read1_sup'] > 0): # only output if support in both strands
                    output_vals = [contig] + list(var) + [v['f_read1_sup'], f_read1_cov, v['r_read1_sup'], r_read1_cov, v['f_sup'], f_cov, v['r_sup'], r_cov, \
                                                          v['f_read1_conc'], v['r_read1_conc'], v['mq'], v['bq'], max_mr]
                    output_file.write('\t'.join(map(str, output_vals)) + '\n')
                del variants[var]

            # flush coverage dictionaries
            to_flush = []
            for frag in frag_coverage:
                if frag[0] + frag[1] < cur_pos: # if frag_start + frag_len is before the current read, we won't need it any more
                    to_flush.append(frag)
                else:
                    break # since frag_coverage is sorted by start position, don't bother checking the remaining fragments (though some could be flushed if they're inside a long fragment)
            for frag in to_flush:
                del frag_coverage[frag]

            ready_to_flush = False

    pbar.close()
    output_file.close()
    tqdm.write(f'finished calling {contig}\n', file=sys.stderr)
    sys.stderr.flush()


# In[ ]:


if __name__ == '__main__':
    chroms = [x.contig for x in aln.get_index_statistics()]
    chroms.sort()

    # start worker processes from pool
    sys.stderr.write('Iterating over chromosomes\n')
    tmp_files = []
    tmp_base = args.tmp + 'tmp_' + args.output.split('/')[-1]
    with multiprocessing.Pool(processes=max(1, args.threads - 1)) as pool:
        lock = multiprocessing.Manager().Lock()
        processes = []
        stats = sorted(aln.get_index_statistics(), key=lambda x: x.contig) # sort contigs by name, this is needed for duplex_caller_output to simultaneously parse multiple variant tsvs
        for i, x in enumerate(stats): # for each chromosome, add a process
            tmp_files.append('{}_{}.tsv'.format(tmp_base, x.contig))
            pool_args = (args.input, x.contig, x.total, aln.get_reference_length(x.contig), tmp_files[-1], i, lock)
            processes.append(pool.apply_async(call_variants, pool_args))
        pool.close()
        pool.join()
    sys.stderr.flush()
    time.sleep(2)
    sys.stderr.write('Checking for errors during runtime\n')
    # check for errors
    for p in processes:
        p.get()


# In[ ]:


# merge and delete temporary files
sys.stderr.write('Merging temporary output files\n')
out_header = '\t'.join(output_columns)
if args.output:
    with open(args.output, 'w') as out_file:
        out_file.write(out_header + '\n')
        for tmp in tmp_files:
            with open(tmp, 'r') as f:
                for l in f:
                    out_file.write(l)
else:
    sys.stdout.write(out_header)
    for tmp in tmp_files:
        with open(tmp, 'r') as f:
            for l in f:
                sys.stdout.write(l) # FIXME this may give a broken pipe error


# In[ ]:


sys.stderr.write('Deleting temporary files\n')
for tmp in tmp_files:
    os.remove(tmp)


# In[ ]:


sys.stderr.write('Complete\n')


# In[ ]:




