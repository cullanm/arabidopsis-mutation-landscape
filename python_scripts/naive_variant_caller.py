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
from os.path import exists
import os
import time
import random
import gtools


# In[ ]:


parser = argparse.ArgumentParser(description='''
Calls variants from a coordinate sorted BAM file without any fancy realignment stuff. Output will treat back-to-back SNVs 
as multiple 1bp variants. Indels of length >1 bp will not be split up. Outputs a tsv with columns:|n
    - chrom: chromosome/contig|n
    - pos: 0-based position of the variant. For insertions, this is the site in the reference which the insertion occurs before|n
    - ref: reference sequence, "*" for insertion|n
    - alt: alt sequence, "*" for deletion|n
    - f_sup: number of forward reads supporting the variant|n
    - f_cov: number of forward reads covering the variant position|n
    - r_sup: number of reverse reads supporting the variant|n
    - r_cov: number of reverse reads covering the variant position|n
    - mapping quality: average mapping quality of supporting reads. If -a is set, ascii converted (qual + 33) phred scaled mapping qualities of supporting reads. Length equals alt depth|n
    - base quality: average base call qualities of supporting bases. If -a is set, ascii converted (qual + 33) phred scaled base call qualities of supporting bases. Length equals alt
        depth * length of variant (e.g. a variant *->TA with two supporting reads where the base call qualities for the two reads
        are ",," and "<<" will have a base quality of ",,<<")
    ''', formatter_class=gtools.EscapeFormatter)
parser.add_argument(dest='input', metavar='INPUT_FILE', type=str, 
                   help='BAM file to call from. Must be coordinate sorted and have an index file')
parser.add_argument(dest='output', metavar='OUTPUT_FILE', type=str,
                   help='output variant file')
parser.add_argument('-@', '--threads', dest='threads', metavar='THREADS_INT', type=int, default=1,
                   help='chromosomes are divided between threads. Cannot use more threads than num chromosomes + 1 (default: 1)')
parser.add_argument('-t', '--tmp', dest='tmp', metavar='TMP_STRING', type=str, default='./',
                   help='prefix to append temporary output files (default: ./)')
parser.add_argument('-q', '--min_mapq', dest='min_mapq', metavar='MIN_MAPQ', type=int, default=1,
                   help='''minimum read mapping quality required for a read to be processed. Reads with lower mapping quality
                    will not count towards any columns (default: 1)''')
parser.add_argument('-s', '--min_support', dest='min_support', metavar='MIN_SUPPORT', type=int, default=1,
                   help='minimum number of supporting reads (forward and reverse) required to output a variant (default: 1)')
parser.add_argument('-p', '--positions', dest='positions', action='store_true', 
                   help='report the position of the variant within all supporting reads as the last column of output (default: do not report positions)')
parser.add_argument('-a', '--all_qualities', dest='all_qualities', action='store_true', 
                   help='report every quality value in the mq and bq column rather than the average')
parser.add_argument('-C', '--chunk_size', dest='chunk_size', metavar='INT', type=int, default=20000,
                   help='''No effect on output. Number of reads to process before flushing buffered information. Increasing
                   may speed up runtime, especially for high-coverage samples. Decreasing may lower memory requirement (default: 20000)''')

try: # run this if in a jupyter notebook
    get_ipython()
    # args = parser.parse_args('-@ 1 -p ../tmp/test.bam ../tmp/test.tsv'.split()) # used for testing
    args = parser.parse_args('-@ 1 -p ../../TEST_OxK_ACR_labelsKi3_chr10%3A100539221-100539618.bam ../../TEST_new_output.tsv'.split()) # used for testing
    # args = parser.parse_args('-@ 1 -a tests/test_naive_variant_caller/test_case.bam tests/test_naive_variant_caller/output.tsv'.split()) # used for testing
    # args = parser.parse_args('../data/align/big_ku80_1_filtered.bam ../tmp/test_naive_variant_caller.tsv -@ 6 --all_qualities -s 2'.split()) # used for testing
except: # run this if in a terminal
    args = parser.parse_args()

if not exists(args.input):
    sys.stderr.write('ERROR: input file {} does not exist\n'.format(args.input))
    exit()

args.output = args.output.strip()
if len(args.output) == 0:
    sys.stderr.write('ERROR: output file is an empty string\n')
    exit()
if args.output[-1] == '/':
    sys.stderr.write('ERROR: output must be a file, not a directory\n')
    exit()

aln = pysam.AlignmentFile(args.input, 'rb') # load BAM for doing checks

# assert index can be opened
try:
    aln.check_index()
except:
    sys.stderr.write('ERROR: failed to open index file {}\n'.format(args.input + '.bai'))
    exit()

# assert output file can be opened
if '/' in args.output and not os.path.isdir(args.output.rsplit('/', 1)[0]):
    sys.stderr.write('ERROR: failed to open output directory for file {}\n'.format(args.output))
    exit()

# assert tmp dir can be opened
if '/' in args.tmp and not os.path.isdir(args.tmp.rsplit('/', 1)[0]):
    sys.stderr.write('ERROR: failed to open tmp directory {}\n'.format(args.tmp.rsplit('/', 1)[0]))
    exit()

# set threads to num chromosomes + 1 if more were allocated
if args.threads > len(aln.get_index_statistics()) + 1:
    args.threads = len(aln.get_index_statistics()) + 1
    sys.stderr.write('WARNING: specified more threads than number of chromosomes + 1. Reducing threads to {}\n'.format(len(aln.get_index_statistics()) + 1))
    
sys.stderr.write('Running naive_variant_caller.py with arguments:\n' + '\n'.join([f'{key}={val}' for key, val in vars(args).items()]) + '\n')


# In[ ]:


# confirm that file is coordinate sorted and contains the MD tag
last = ('', 0)
for i, read in enumerate(aln.fetch()):
    nex = (read.reference_name, read.reference_start)
    if nex < last:
        sys.stderr.write('ERROR: SAM/BAM file is not coordinate sorted\n')
        exit()
    last = nex
    
    try:
        read.get_tag('MD')
    except KeyError:
        sys.stderr.write('ERROR: SAM/BAM file does not contain the "MD" tag\n')
        exit()
    
    if i > 1000:
        break


# ### Notes:
# - code requires coordinate sorted BAM
# - code requires the SAM to contain the MD tag to function (see https://samtools.github.io/hts-specs/SAMtags.pdf). Tag follows the regex `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`, where a number indicates a match with the reference, a letter indicates the reference base at a mismatch, and a deletion is specified as a carat followed by one or more letters indicating the aligned reference sequence. This is how Bowtie2 produces the MD field. Insertions are not represented in the MD tag.
# - positions are 0 based in varaibles "vcf" and "variants"

# In[ ]:


draw_pbars = True # whether to draw tqdm progress bars


# In[ ]:


output_cols = '\t'.join('chrom pos ref alt f_sup f_cov r_sup r_cov mq bq read_pos'.split())


# In[ ]:


# calls variants for a single contig. reads must be an iterator over a contig from AlignmentFile.fetch()
def call_variants(input_filepath, contig, total_reads, contig_len, output_filepath, process_num, lock):
    aln = pysam.AlignmentFile(input_filepath, 'rb')

    # key is (pos (0 based), ref, alt)
    variants = defaultdict(lambda: {'f_sup': 0, # supporting forward reads
                                    'f_cov': 0, # coverage of forward reads
                                    'r_sup': 0, # supporting reverse reads
                                    'r_cov': 0, # coverage of reverse reads
                                    'mq': '', # mapping qualities of supporting reads
                                    'bq': '', # base qualities of supporting reads
                                    'positions': []}) # position of variant in supporting reads
    f_coverage = np.zeros(contig_len, dtype=np.intc) # forward strand coverage
    r_coverage = np.zeros(contig_len, dtype=np.intc) # reverse strand coverage
    num_read = 0
    cur_pos = -1
    ready_to_flush = False
    output_file = open(output_filepath, 'w') # temporary output file
    
    # initialize progress bar
    if draw_pbars:
        with lock:
            pbar = tqdm(total=total_reads, desc=contig, position=process_num, leave=True, miniters=100000, maxinterval=9999)
    
    for read in aln.fetch(contig=contig):
        num_read += 1
        
        # skip pairs that don't pass the mapq filter
        if read.mapq < args.min_mapq:
            continue

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
        
        print(mds)
        
        # cigar strings include insertions, soft, and hard clipping, info not present in the MD tag
        # get the length of starting softclips and the number of clipping operations at the start
        softclip_start_len = 0
        clip_cigs = 0
        for cig in read.cigartuples:
            if cig[0] == 5:
                clip_cigs += 1
            elif cig[0] == 4:
                softclip_start_len += cig[1]
                clip_cigs += 1
            else:
                break
                
        softclip_end_len = 0
        for cig in read.cigartuples[-1:-1:-1]:
            if cig[0] == 4:
                softclip_end_len += cig[1]
            elif cig[0] != 5:
                break
        
        cur_cig = 0 # keep track of the sum of cigar element lengths
        # find the insertions in the cigar, add them to mds as "*x", where x in the length of insertion
        for cig in read.cigartuples[clip_cigs:]: # ignore the clipping cigars
            if cig[0] == 1: # if insertion (I)
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
            elif cig[0] != 0 and cig[0] != 2 and cig[0] != 4 and cig[0] != 5: # if cigar operation isn't M, D, I, S, nor H
                with lock:
                    sys.stderr.write('WARNING: treating unexpected cigar operation "{}" as match in read:\n{}\n'.format(cig[0], read))
            cur_cig += cig[1] # note: clipping shouldn't be added here, but it doesn't matter since all remaining clipping is at the end
        
        # remove spacer zeros (e.g. ['A', '0', 'T'] -> ['A', 'T'])
        mds = [md for md in mds if md != '0']
                
        # example of how mds should now look: ['A', '3', 'C', 'G', '^TA', '9', '*2', '2']
        
        has_bqs = read.query_qualities is not None # query_qualities can be none if the BQ field is a sigle "*"
        
        # find and add all variants using md
        cur_md = 0
        total_insertion_len = 0
        for md in mds:
            if md.isnumeric(): # skip over matches
                cur_md += int(md)
            else:
                read_pos = cur_md
                if md[0] == '^': # if deletion
                    key = (read.reference_start + cur_md - total_insertion_len, md[1:], '*')
                    qualities = ''# no base call qualities to add
                    # don't increment cur_md, as deletions don't count towards the read length
                    total_insertion_len -= len(md) - 1
                elif md[0] == '*': # if insertion
                    key = (read.reference_start + cur_md - total_insertion_len, '*', read.query_sequence[softclip_start_len + cur_md:softclip_start_len + cur_md + int(md[1:])])
                    qualities = ''
                    if has_bqs: 
                        for qual in read.query_qualities[cur_md:cur_md+int(md[1:])]:
                            qualities += chr(qual + 33)
                    cur_md += int(md[1:])
                    total_insertion_len += int(md[1:])
                else: # if mismatch
                    key = (read.reference_start + cur_md - total_insertion_len, md, read.query_sequence[softclip_start_len + cur_md])
                    qualities = chr(read.query_qualities[cur_md] + 33) if has_bqs else ''
                    cur_md += 1
                
                if read.is_reverse:
                    variants[key]['r_sup'] += 1
                    variants[key]['positions'].append(read.query_length - read_pos - 1 + softclip_end_len)
                else:
                    variants[key]['f_sup'] += 1
                    variants[key]['positions'].append(read_pos + softclip_start_len)
                variants[key]['mq'] += chr(read.mapping_quality + 33)
                variants[key]['bq'] += qualities
    
        # add to read coverage array
        if read.is_reverse:
            r_coverage[read.reference_start:read.reference_start + read.reference_length] += 1
        else:
            f_coverage[read.reference_start:read.reference_start + read.reference_length] += 1
        
        # periodically flush variants and coverage info out of the dictionaries
        if num_read % args.chunk_size == 0:
            if draw_pbars:
                with lock:
                    pbar.update(args.chunk_size)
            ready_to_flush = True
        
        # periodically flush variants out of the dictionary
        if (ready_to_flush and read.reference_start > cur_pos) or num_read == total_reads:
            # variants = defaultdict(variants.default_factory, {key: val for key, val in sorted(variants.items(), key=lambda ele: ele[0])}) # sort by position
            cur_pos = read.reference_start
            
            to_flush = []
            for var in variants:
                if var[0] >= cur_pos and num_read < total_reads:
                    continue
                to_flush.append(var)
            for var in to_flush:
                if variants[var]['f_sup'] + variants[var]['r_sup'] >= args.min_support:
                    variants[var]['f_cov'] = f_coverage[var[0]]
                    variants[var]['r_cov'] = r_coverage[var[0]]
                    v = variants[var]
                    
                    if not args.all_qualities:
                        v['mq'] = sum([ord(x) - 33 for x in v['mq']]) / len(v['mq'])
                        if var[2] == '*' or v['bq'] == '':
                            v['bq'] = 42
                        else:
                            v['bq'] = sum([ord(x) - 33 for x in v['bq']]) / len(v['bq'])
                    
                    out_line = '\t'.join(map(str, [contig] + list(var) + [v['f_sup'], v['f_cov'], v['r_sup'], v['r_cov'], v['mq'], v['bq']]))
                    out_line += '\t' + ','.join(map(str, v['positions'])) if args.positions else ''
                    output_file.write(out_line + '\n')
                del variants[var]
            if draw_pbars:
                with lock:
                    pbar.update(args.chunk_size)
            ready_to_flush = False

    output_file.close()


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

with open(args.output, 'w') as out_file:
    out_file.write(output_cols + '\n')
    for tmp in tmp_files:
        with open(tmp, 'r') as f:
            for l in f:
                out_file.write(l)


# In[ ]:


sys.stderr.write('Deleting temporary files\n')
for tmp in tmp_files:
    os.remove(tmp)


# In[ ]:


sys.stderr.write('Complete\n')


# In[ ]:





# In[ ]:




