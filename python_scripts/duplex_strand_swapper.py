#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pysam
from tqdm import tqdm
import random
import subprocess
from multiprocessing import Process
from collections import defaultdict
import sys
import os
import subprocess
import multiprocessing
import itertools
import argparse
import gtools


# In[ ]:


parser = argparse.ArgumentParser(description='''
Makes a randomly strand-swapped NanoSeq library from several input NanoSeq libraries. For each fragment start and end position in the input libraries, one
set of top strand PCR duplicates (F1R2) and one set of bottom strand PCR duplicates (F2R1) is randomly selected from two different libraries. If this isn't 
possible, nothing is output for that fragment. The UMI tag will be overwritten as {index of top strand BAM}:{top strand UMI}_{index of bot strand BAM}:{bot strand UMI}.|n
|n
The point of this script is to make all the true mutations uncallable, while keeping the false positive rate the same. Because FPs appear in strands 
independently, it shouldn't matter how we rearrange the strands, the likelihood of pairing two strands with the same FP remains the same. |n
|n
Example:|n 
For one fragment start and end position, three libraries have 2 fragments with reads from both the top and bottom strand. Out of the 6 fragments, we randomly 
select sample_2 umi_2 to contribute the top strand reads. Out of the 4 fragments not in sample_2, we randomly select sample_0 umi_2 to contribute the bottom 
strand reads. One fragment will be output with the reads of sample_2 umi_2's top strand and sample_0 umi_2's bottom strand, and the umi will be 
"2:umi_2_0:umi_2".
    ''', formatter_class=gtools.EscapeFormatter)
parser.add_argument('-i', '--input', required=True, dest='input', metavar='FILES', type=str, nargs='+',
                   help='input NanoSeq BAM files. Must be coordinate sorted and indexed')
parser.add_argument('-o', '--output', required=True, dest='output', metavar='FILE', type=str,
                   help='output strand-swapped BAM file')
parser.add_argument('-u', '--umi', dest='umi', metavar='TAG', type=str, default=None,
                   help='name of tag which contains the fragment umi')
parser.add_argument('-t', '--tmp', dest='tmp_prefix', metavar='PREFIX', type=str, default='tmp/',
                   help='prefix to add to temporary files. Will make any directories which do not exist (default: tmp/)')
parser.add_argument('-@', '--threads', dest='threads', metavar='INT', type=int, default=1,
                   help='number of threads. Each chromosome can be processed by a separate thread. (default: 1)')

try: # run this if in a jupyter notebook
    get_ipython()
    print('Determined code is running in Jupyter')
    args = parser.parse_args('-@ 8 --umi RG -o tmp/swapper_test.bam -i data/align/big_Col-0-1_merged.bam'.split()) # used for testing
except: # run this if in a terminal
    args = parser.parse_args()

sys.stderr.write('Running duplex_strand_swapper.py with arguments:\n' + '\n'.join([f'{key}={val}' for key, val in vars(args).items()]) + '\n')
if args.tmp_prefix and '/' in args.tmp_prefix:
    os.makedirs(os.path.dirname(args.tmp_prefix), exist_ok=True)

if len(args.input) == 1:
    sys.stderr.write('ERROR: must specificy more than one input file\n')
    sys.exit(1)


# In[ ]:


update_interval = 1000000


# In[ ]:


def strand_swap(output_file, chrom):
    stream_in = gtools.BamStreamer(args.input, chrom) # iterator over reads in multiple files, returns (file index, read) tuples
    aln_out = pysam.AlignmentFile(output_file, 'wb', template=stream_in.aln_files[0])
    cur_source, cur_read = next(stream_in)
    reverse_to_output = dict() # whenever a forward read is output, add an item of query name:new RG tag
    done = False
    next_update = update_interval
    while not done:
        cur_pos = cur_read.reference_start

        # get all the reads with the same starting position
        reads = defaultdict(lambda: []) # stores the forward reads, key=tlen, value=[(file source, read)]
        reverse_reads = [] # stores the reverse reads
        while cur_read.reference_start == cur_pos:
            if cur_read.is_reverse:
                reverse_reads.append(cur_read)
            else:
                reads[abs(cur_read.template_length)].append((cur_source, cur_read))
            try:
                cur_source, cur_read = next(stream_in)
            except StopIteration:
                done = True
                break


        for tlen in reads:
            top_sources = list(set([(source, read.get_tag(args.umi)) for source, read in reads[tlen] if read.is_read1]))
            bot_sources = list(set([(source, read.get_tag(args.umi)) for source, read in reads[tlen] if read.is_read2]))

            if len(top_sources) == 0:
                continue

            top_sample = random.choice(top_sources)
            bot_sources = [(source, umi) for source, umi in bot_sources if source != top_sample[0]]

            if len(bot_sources) == 0:
                continue

            bot_sample = random.choice(bot_sources)

            for source, read in reads[tlen]:
                if (read.is_read1 and source == top_sample[0] and read.get_tag(args.umi) == top_sample[1]) or \
                   (read.is_read2 and source == bot_sample[0] and read.get_tag(args.umi) == bot_sample[1]):
                    read.set_tag(args.umi, f'{top_sample[0]}:{top_sample[1]}_{bot_sample[0]}:{bot_sample[1]}')
                    aln_out.write(read)
                    reverse_to_output[read.query_name] = f'{top_sample[0]}:{top_sample[1]}_{bot_sample[0]}:{bot_sample[1]}'

        # output any reverse reads 
        for read in reverse_reads:
            if read.query_name in reverse_to_output:
                read.set_tag(args.umi, reverse_to_output[read.query_name])
                aln_out.write(read)
                del reverse_to_output[read.query_name]

        if cur_pos >= next_update:
            sys.stderr.write(f'passed {chrom}:{next_update}\n')
            next_update += update_interval

    stream_in.close()
    aln_out.close()


# In[ ]:


a = pysam.AlignmentFile(args.input[0], 'rb')
chroms = a.references
a.close()


# In[ ]:


if __name__ == '__main__':
    # start worker processes from pool
    pool = multiprocessing.Pool(processes=args.threads)
    processes = []
    for chrom in chroms:
        processes.append(pool.apply_async(strand_swap, (f'{args.tmp_prefix}{args.output.replace("/", ".")}_{chrom}.bam', chrom)))
    pool.close()
    pool.join()

    # check for errors
    for p in processes:
        p.get()


# In[ ]:


to_cat = ' '.join([f'{args.tmp_prefix}{args.output.replace("/", ".")}_{chrom}.bam' for chrom in chroms])

p = subprocess.run(f'samtools cat -o {args.output} {to_cat}', shell=True, capture_output=True)
if p.returncode != 0:
    sys.stderr.write(str(p.stderr))
    exit()

p = subprocess.run(f'rm {to_cat}', shell=True, capture_output=True)
if p.returncode != 0:
    sys.stderr.write(str(p.stderr))
    exit()


# In[ ]:


sys.stdout.write(f'completed duplex_strand_swapper with output {args.output}\n')


# In[ ]:




