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
Calls variants from a coordinate sorted BAM file without any fancy realignment stuff. Ignores reads that didn't map as a concordant pair. Output will treat back-to-back SNVs 
as multiple 1bp variants. Indels of length >1 bp will appear as a single row. Takes pcr/optical duplicates and read1-read2 overlap into consideration when determining
if a read is in support and doesn't double count duplicates. For each set of reads with the same fragment start and end, they only contribute to support if variants are supported by > x% of the covering 
reads in the set (includes both forward and reverse reads). If a set of reads passes this filter, it contributes 1 to forward/reverse support (or both if in the 
read1-read2 overlap). This means read pairs with no duplicates still must pass this filter within the read1-read2 overlap. Outputs a tsv with columns:|n
    - chrom: chromosome/contig|n
    - pos: 0-based position of the variant. For insertions, this is the site in the reference which the insertion occurs before|n
    - ref: reference sequence, "*" for insertion|n
    - alt: alt sequence, "*" for deletion|n
    - f_sup: number of forward reads supporting the variant|n
    - f_cov: number of forward reads covering the variant position|n
    - r_sup: number of reverse reads supporting the variant|n
    - r_cov: number of reverse reads covering the variant position|n
    - mapping quality: ascii converted (qual + 33) phred scaled mapping qualities of supporting reads. Length equals alt depth. Selects the median value for sets of duplicates|n
    - base quality: ascii converted (qual + 33) phred scaled base call qualities of bases in supporting reads. Length equals alt
        depth * length of variant (e.g. a variant *->TA with two supporting reads where the base call qualities for the two reads
        are ",," and "<<" will have a base quality of ",,<<"). Selects the median value for sets of duplicates
    ''', formatter_class=gtools.EscapeFormatter)
parser.add_argument(dest='input', metavar='FILE', type=str, 
                   help='BAM file to call from. Must be coordinate sorted and have an index file')
parser.add_argument(dest='output', metavar='FILE', type=str,
                   help='output variant file')
parser.add_argument('-@', '--threads', dest='threads', metavar='INT', type=int, default=1,
                   help='chromosomes are divided between threads. Cannot use more threads than num chromosomes + 1 (default: 1)')
parser.add_argument('-t', '--tmp', dest='tmp', metavar='PREFIX', type=str, default='./',
                   help='prefix to append temporary output files (default: ./)')
parser.add_argument('-q', '--min_mq', dest='min_mq', metavar='INT', type=int, default=1,
                   help='''minimum read mapping quality required for a read to be processed. Reads with lower mapping quality
                    will not count towards any columns (default: 1)''')
parser.add_argument('-d', '--duplicate_ratio', dest='duplicate_ratio', metavar='FLOAT', type=float, default=0.79,
                   help='''Cutoff for whether a group of duplicates supports a variant. If the number of supporting reads from
                   all duplicates / the number of covering reads from all duplicates (both forward and reverse) is greater than 
                   DUPLICATE_RATIO, the group of duplicates will add 1 to the variant's support (1 to forward and reverse 
                   support if in the read1-read2 overlap) (default: 0.79)''')
parser.add_argument('-p', '--duplicate_support', dest='duplicate_sup', metavar='INT', type=int, default=1,
                   help='''Minimum number of duplicates required to contibute to variant support and coverage (default: 1)''')
parser.add_argument('-m', '--min_support', dest='min_support', metavar='INT', type=int, default=1,
                   help='only output variants with f_sup + r_sup > MIN_SUPPORT (default: 1)')
parser.add_argument('-s', '--min_strand_support', dest='min_strand_support', metavar='INT', type=int, default=0,
                   help='only output variants with f_sup > MIN_STRAND_SUPPORT and r_sup > MIN_STRAND_SUPPORT (default: 0)')
parser.add_argument('-u', '--umi', dest='umi', metavar='TAG', type=str, default=None,
                   help='Consider the UMI in the provided BAM tag when grouping pre-PCR fragments (default: no UMI)')
parser.add_argument('-e', '--end_dist', dest='end_dist', action='store_true',
                   help='Report the distance from fragment end of each supporting read as an extra column (default: not reported)')
parser.add_argument('-C', '--chunk_size', dest='chunk_size', metavar='INT', type=int, default=50000,
                   help='''No effect on output. Number of reads to process before flushing buffered information. Increasing
                   may speed up runtime, especially for high-coverage samples. Decreasing lowers memory requirement (default: 10000)''')

try: # run this if in a jupyter notebook
    get_ipython()
    if os.getcwd()[:8] != '/scratch': # switch to the scratch directory where all the data files are
        os.chdir(f'/scratch/cam02551/{os.getcwd().split("/")[-2]}')
    # args = parser.parse_args('tests/test_dup_informed_caller/dup_informed_caller_test_case.bam tests/test_dup_informed_caller/test.out -@ 1 -C 1 -p 4'.split()) # used for testing
    args = parser.parse_args('data/align/chandler_KVKCS001D_0_filtered.bam tmp/dup_inf_testing.tsv --umi RX --duplicate_support 2 -@ 4'.split()) # used for testing
except: # run this if in a terminal
    args = parser.parse_args()

aln = pysam.AlignmentFile(args.input, 'rb') # load BAM for doing checks

# assert index can be opened
try:
    aln.check_index()
except:
    sys.stderr.write('ERROR: failed to open index file {}\n'.format(args.input + '.bai'))
    exit()

# assert output file can be opened
if args.output and '/' in args.output and not os.path.isdir(args.output.rsplit('/', 1)[0]):
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

sys.stderr.write('Running dup_informed_caller.py with arguments:\n' + '\n'.join([f'{key}={val}' for key, val in vars(args).items()]) + '\n')


# In[ ]:


# confirm that file is coordinate sorted and contains the MD tag and specified UMI
last = ('', 0)
# max_read_length = 0
for i, read in enumerate(aln.fetch()):
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
    
    # max_read_length = max(max_read_length, read.query_length)
    if i > 1000:
        break
# sys.stderr.write(f'inferring a read length of {max_read_length}bp\n')


# ### Notes:
# - code requires coordinate sorted BAM
# - code requires the BAM to contain the MD tag to function (see https://samtools.github.io/hts-specs/SAMtags.pdf). Tag follows the regex `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`, where a number indicates a match with the reference, a letter indicates the reference base at a mismatch, and a deletion is specified as a carat followed by one or more letters indicating the aligned reference sequence. This is how Bowtie2 produces the MD field. Insertions are not represented in the MD tag.
# - all position variables and output should be zero-based

# In[ ]:


draw_pbars = True # whether to draw tqdm progress bars


# In[ ]:


# these will be the column names ofthe final output, the temporary outputs of call_variants must match this order
output_columns = 'chrom pos ref alt f_sup f_cov r_sup r_cov mq bq'.split()


# In[ ]:


# calls variants for a single contig
def call_variants(input_filepath, contig, total_reads, contig_len, output_filepath, process_num, lock):
    aln = pysam.AlignmentFile(input_filepath, 'rb')

    # key is (pos (0 based), ref, alt, frag_start, frag_length, umi)
    variants = defaultdict(lambda: {'f_sup': 0, # supporting reads on + strand
                                    'r_sup': 0, # supporting reads on - strand
                                    'mq': '', # mapping qualities of supporting reads
                                    'bq': ''}) # base qualities of supporting reads
    frag_coverage = dict() # key is (frag_start, frag_length, strand, umi); value is array of read coverage over the fragment
    uncounted_frags = [] # frags in frag_coverage which haven't yet contributed to chrom_coverage
    chrom_coverage = np.zeros((2, contig_len), dtype=np.intc) # forward [0] and reverse [1] coverage of contig (duplicates don't count)
    num_read = 0
    cur_pos = -1
    ready_to_flush = False
    output_file = open(output_filepath, 'w') # temporary output file
    
    # initialize progress bar
    if draw_pbars:
        with lock:
            pbar = tqdm(total=total_reads, desc=contig, position=process_num, leave=True, mininterval=10, maxinterval=9999)
    
    for read in aln.fetch(contig=contig):
        num_read += 1
        
        # skip pairs that don't align as expected or don't pass the mq filter. this might cause incorrect coverage values if two mates have different mqs
        if read.is_unmapped or read.mate_is_unmapped or not (read.flag & 2) or read.mapping_quality < args.min_mq:
            continue
            
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
        
        # get the UMI tag
        if args.umi is None:
            umi = ''
        else:
            umi = read.get_tag(args.umi) 
        
        # find and add all variants using md
        cur_md = 0 # current position in the read
        total_insertion_len = 0 # aligned read is
        for md in mds:
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
                       frag_start, frag_len, umi)
                    
                # add variants to dictionary
                if read.is_reverse:
                    variants[key]['r_sup'] += 1
                else:
                    variants[key]['f_sup'] += 1
                variants[key]['mq'] += chr(read.mapping_quality + 33)
                variants[key]['bq'] += qualities
        
        # add coverage of fragment to frag_coverage (and to chrom_coverage if a duplicate wasn't counted yet)
        frag = (frag_start, frag_len, umi)
        
        # if this is the first read seen of a fragment, add it to the frag_coverage dict
        if frag not in frag_coverage: 
            frag_coverage[frag] = np.zeros((2, frag[1]), dtype=np.short)
            uncounted_frags.append(frag)
        
        # add the coverage of the read and it's mate to frag_coverage
        if not read.is_reverse: # only run for forward reads, since we find these first and can infer reverse read coverage using read.next_reference_start
            frag_coverage[frag][0, :read.reference_length] += 1
            frag_coverage[frag][1, read.next_reference_start - read.reference_start:] += 1
            
        # periodically flush variants and coverage info out of the dictionaries
        if num_read % args.chunk_size == 0:
            if draw_pbars:
                with lock:
                    pbar.update(args.chunk_size)
            ready_to_flush = True

        if (ready_to_flush and read.reference_start > cur_pos) or num_read == total_reads:
            # variants = defaultdict(variants.default_factory, {key: val for key, val in sorted(variants.items(), key=lambda ele: ele[0])}) # sort by position, not necessary
            cur_pos = read.reference_start
            
            # compute chromosome level coverage up to cur_pos
            last_counted_frag = -1
            for i, f in enumerate(uncounted_frags):
                # if we reach the current read position, skip, as we may not have found all the reads for a fragment yet
                if f[0] >= cur_pos: 
                    break
                total_frag_cov = np.sum(frag_coverage[f], axis=0) # summed forward and reverse coverage over the fragment
                high_cov_mask = total_frag_cov >= args.duplicate_sup # positions in the fragment with enough coverage to pass the duplicate_sup requirement, only these positions should contribute to coverage
                high_cov_forward = np.logical_and(high_cov_mask, frag_coverage[f][0]) # positions with enough coverage and >0 coverage from forward reads
                high_cov_reverse = np.logical_and(high_cov_mask, frag_coverage[f][1]) # same but for reverse reads

                # add positions with sufficient coverage to chrom_coverage
                chrom_coverage[0, f[0]:f[0] + f[1]] += high_cov_forward
                chrom_coverage[1, f[0]:f[0] + f[1]] += high_cov_reverse
                
                last_counted_frag = i
            uncounted_frags = uncounted_frags[last_counted_frag + 1:] # remove all counted fragments
            
            # group variants from multiple fragments
            to_flush = dict() # key is (pos, ref, alt) value is [f_sup, r_sup]
            to_del = set()
            for var in variants: # (pos, ref, alt, frag_start, frag_len)
                if var[0] >= cur_pos and num_read < total_reads:
                    continue
                v = variants[var]
                f_cov = frag_coverage[(var[3], var[4], var[5])][0, var[0] - var[3]]
                r_cov = frag_coverage[(var[3], var[4], var[5])][1, var[0] - var[3]]
                if (v['f_sup'] + v['r_sup']) / (f_cov + r_cov) > args.duplicate_ratio and v['f_sup'] + v['r_sup'] >= args.duplicate_sup: # if it passes the duplicate ratio and support filter
                    glob_var = (var[0], var[1], var[2])
                    if glob_var not in to_flush:
                        to_flush[glob_var] = [0, 0, '', '']
                    to_flush[glob_var][0] += 1 if f_cov > 0 else 0
                    to_flush[glob_var][1] += 1 if r_cov > 0 else 0
                    to_flush[glob_var][2] += sorted(v['mq'])[len(v['mq']) // 2] # take median mq
                    
                    # take median bq(s)
                    if var[1] == '*': # if insertion
                        for i in range(len(var[2])): # get the median bq of each position (e.g. for a bq of 123456, output median(1,3,5) median(2,4,6))
                            to_flush[glob_var][3] += sorted(v['bq'][i::len(var[2])])[len(v['mq']) // 2]
                    elif var[2] == '*': # if deletion
                        pass
                    else: # if SNV
                        to_flush[glob_var][3] += sorted(v['bq'])[len(v['bq']) // 2]
                to_del.add(var)

            for var in to_del:
                del variants[var]

            # output the variants
            for var in sorted(to_flush): # (pos, ref, alt)
                f_cov = chrom_coverage[0, var[0]]
                r_cov = chrom_coverage[1, var[0]]
                sup = to_flush[var]
                
                # check if enough support
                if (sup[0] >= args.min_strand_support) and (sup[1] >= args.min_strand_support) and (sup[0] + sup[1] >= args.min_support):
                    output_vals = [contig] + list(var) + [sup[0], f_cov, sup[1], r_cov, sup[2], sup[3]]
                    output_file.write('\t'.join(map(str, output_vals)) + '\n')
            # flush coverage dictionaries
            to_flush = []
            for frag in frag_coverage:
                if frag[0] + frag[1] < cur_pos: # if frag_start + frag_len is before the current read, we won't need it any more
                    to_flush.append(frag)
                    if frag in uncounted_frags:
                        print('ruh roh')
            for frag in to_flush:
                del frag_coverage[frag]
            ready_to_flush = False
    
    pbar.close()
    output_file.close()
    tqdm.write(f'finished calling {contig}\n', file=sys.stderr)
    sys.stderr.flush()


# In[ ]:


if __name__ == '__main__':
    # start worker processes from pool
    sys.stderr.write('Running a subprocess for each chromosome\n')
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




