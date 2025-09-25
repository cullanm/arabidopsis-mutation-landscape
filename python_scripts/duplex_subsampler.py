#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import os
import gtools

import pysam
from tqdm import tqdm
import random
import subprocess
from multiprocessing import Process
from collections import defaultdict
import sys
import subprocess
import multiprocessing
import argparse


# In[ ]:


parser = argparse.ArgumentParser(description='''
Subsamples fragments/reads from a NanoSeq library BAM.
    ''')
parser.add_argument(dest='input', metavar='FILE', type=str,
                   help='input NanoSeq BAM file. Must be coordinate sorted and indexed')
parser.add_argument('-o', '--output', required=True, dest='output', metavar='FILE', type=str,
                   help='subsampled BAM file')
parser.add_argument('-s', '--required_strand_coverage', dest='strand_req', metavar='INT', type=int, default=0,
                   help='number of reads in each strand required to output the fragment (default: 0)')
parser.add_argument('-t', '--required_total_coverage', dest='total_req', metavar='INT', type=int, default=0,
                   help='number of total reads in the fragment required to output (default: 0)')
parser.add_argument('-c', '--cap', dest='cap_at_req', action='store_true',
                   help='caps the number of reads output for a fragment at the minimum required to pass the requirements \
                   reads to be output are selected at random without replacement (default: output all reads for passing fragments)')
parser.add_argument('-r', '--random_frac', dest='random_frac', metavar='FLOAT', type=float, default=1.001,
                   help='fraction of fragments passing requirements to output (default: 1 (output all passing fragments))')
parser.add_argument('-u', '--umi', dest='umi', metavar='TAG', type=str, default=None,
                   help='Consider the UMI in the provided BAM tag when grouping fragments as PCR dups (default: no UMI)')
parser.add_argument('--tmp', dest='tmp_prefix', metavar='PREFIX', type=str, default='tmp/',
                   help='prefix to add to temporary files. Will make any directories which do not exist (default: tmp/)')
parser.add_argument('-@', '--threads', dest='threads', metavar='INT', type=int, default=1,
                   help='number of threads. Each chromosome can be processed by a separate thread. (default: 1)')

try: # run this if in a jupyter notebook
    get_ipython()
    print('Determined code is running in Jupyter')
    args = parser.parse_args('-@ 8 --umi RG -s 4 -t 10 --cap -o tmp/subsampler_test.bam tmp/smol.bam'.split()) # used for testing
    if os.getcwd()[:8] != '/scratch': # switch to the scratch directory where all the data files are
        os.chdir(f'/scratch/cam02551/{os.getcwd().split("/")[-2]}')
except: # run this if in a terminal
    args = parser.parse_args()

os.makedirs(os.path.dirname(args.tmp_prefix), exist_ok=True)


# In[ ]:


buffer_size = 200000 # amount of reads to process before each flush
read_len = 150


# In[ ]:


# subsamples one chromosome of a bam file to only output fragments with at least a certain amount of read coverage
# file_in and file_out are bam files
# f_req and r_req are the required number of read pairs from the original forward and reverse strands
# to output only f_req forward and r_req reverse read pairs, set cap_at_req to True
# if cap_at_req is True, reads will be selected at random WITHOUT REPLACEMENT for output
def subsample_frags(file_in, file_out, chrom):
    fragments = defaultdict(lambda: [0, 0, [], []]) # index as fragments[(frag_start, tlen, umi)] = [f_read1_supp, r_read1_supp, [f read names], [r read names]]
    aln_in = pysam.AlignmentFile(file_in, 'rb')
    aln_out = pysam.AlignmentFile(file_out, 'wb', template=aln_in)
    buffer = [] # stores reads in the order read from aln_in, these are eventually written out
    num_written = 0 # reads written
    num_frags = 0 # fragments written
    buffer_log = [] # history of buffer sizes before and after each flush
    to_output = dict() # read names to output from buffer, value is number of copies to output
    to_discard = set() # read names to discard from buffer
    num_read = 0
    for x in aln_in.get_index_statistics():
        if x.contig == chrom:
            total_reads = x.total
    
    for read in tqdm(aln_in.fetch(contig=chrom), mininterval=10, total=total_reads, desc=chrom):
        num_read += 1
        
        # periodically flush out reads
        if num_read % buffer_size == 0 or num_read == total_reads:
            buffer_log.append(len(buffer))
#             print(len(buffer))
            
            new_fragments = defaultdict(lambda: [0, 0, [], []]) # will later overwrite fragments
            for x in fragments:
                
                # if the fragment start is past the reader position, don't consider outputting it yet, as more PCR dups may be found
                if num_read != total_reads and x[0] >= frag_start:
                    new_fragments[x] = fragments[x]
                    
                # if there are enough reads in the fragment to meet the requirements, and the fragment is ranomly selected to output (-r)
                elif fragments[x][0] >= args.strand_req // 2 and fragments[x][1] >= args.strand_req // 2 and fragments[x][0] + fragments[x][1] >= args.total_req // 2 and random.random() < args.random_frac:
                    num_frags += 1
                    if args.cap_at_req: # if subsampling PCR duplicates
                        # choose read names from each original strand without replacement
                        
                        chosen_names = random.sample(fragments[x][2], k=args.strand_req // 2)
                        chosen_names += random.sample(fragments[x][3], k=args.strand_req // 2)
                        if args.total_req > args.strand_req * 2:
                            remaining_names = [name for name in fragments[x][2] + fragments[x][3] if name not in chosen_names]
                            chosen_names += random.sample(remaining_names, k=(args.total_req - 2 * args.strand_req) // 2)
                            
                        for name in chosen_names: # add names to to_output
                            if name in to_output:
                                to_output[name] += 1
                            else:
                                to_output[name] = 1
                        for name in fragments[x][2] + fragments[x][3]: # add unchosen names to to_discard
                            if name not in chosen_names:
                                to_discard.add(name)
                    else: # if not subsampling PCR duplicates, add all reads to to_output
                        to_output.update({name:1 for name in fragments[x][2]})
                        to_output.update({name:1 for name in fragments[x][3]})
                else: # if there aren't enough covering fragments, add reads to to_discard
                    to_discard.update(fragments[x][2])
                    to_discard.update(fragments[x][3])
        
            new_buffer = [] # will later overwrite buffer
            done = False # when to stop outputting reads
            for x in buffer: # for each read in buffer
                if done: # do nothing if done outputting
                    new_buffer.append(x)
                elif x.query_name in to_output: # if read should be output
                    num_written += 1
                    aln_out.write(x)
                    num_copies = to_output[x.query_name] # how many copies of the read to write
                    basename = x.query_name
                    for i in range(1, num_copies): # write extra copies with ":n" appended to the read name
                        x.query_name = basename + f':{i}'
                        num_written += 1
                        aln_out.write(x)
                    if x.is_reverse: # remove from to_output once we output the reverse read
                        del to_output[basename]
                elif x.query_name in to_discard: # if read should be discarded
                    if x.is_reverse: # remove from to_output once we output the reverse read
                        to_discard.remove(x.query_name)
                else: # if don't know whether the read should be output yet, stop outputting to preserve sort order
                    done = True
                    new_buffer.append(x)
            
            fragments = new_fragments
            buffer = new_buffer
        
            buffer_log.append(len(buffer))

        
        # thow out read pairs that don't align as expected
        if read.is_unmapped or read.mate_is_unmapped or not (read.flag & 2):
            continue
        
        buffer.append(read)
        
        if read.is_reverse: # ignore reverse reads so we don't double count each fragment
            continue
        
        # find the fragment start and end. Used to find PCR duplicates
        frag_start = min(read.reference_start, read.next_reference_start)
        tlen = abs(read.template_length)
        umi = read.get_tag(args.umi)
        
        read_strand = '-' if read.is_read2 else '+'
        
        # add to the fragments read pair count and read name list
        if read_strand == '+':
            fragments[(frag_start, tlen, umi)][0] += 1
            fragments[(frag_start, tlen, umi)][2].append(read.query_name)
        else:
            fragments[(frag_start, tlen, umi)][1] += 1
            fragments[(frag_start, tlen, umi)][3].append(read.query_name)
        
        # FIXME this loop freezes on some chromosomes (e.g. Chr3 52%, Chr4 64%, Chr3 52%)
    aln_in.close()
    aln_out.close()
    sys.stderr.write(f'completed chromosome {chrom}, wrote {num_written} reads and {num_frags} fragments,     max buffer size of {max(buffer_log)} reads at {buffer_log.index(max(buffer_log)) / len(buffer_log) * 100:.1f}%\n')


# # Make a file containing all reads from all high coverage fragments

# In[ ]:


aln = pysam.AlignmentFile(args.input)
chroms = aln.references
aln.close()


# In[ ]:


sys.stderr.write('subsampling each chromosome\n')
if __name__ == '__main__':
    # start worker processes from pool
    pool = multiprocessing.Pool(processes=args.threads)
    processes = []
    for chrom in chroms:
        processes.append(pool.apply_async(subsample_frags, (args.input, f'{args.tmp_prefix}{args.output.replace("/", ".")}_{chrom}.bam', chrom)))
    pool.close()
    pool.join()

    # check for errors
    for p in processes:
        p.get()


# In[ ]:


sys.stderr.write('concatenating subsampled chromosome BAMs\n')

to_cat = ' '.join([f'{args.tmp_prefix}{args.output.replace("/", ".")}_{chrom}.bam' for chrom in chroms])

p = subprocess.run(f'samtools cat -@ {args.threads} -o {args.output} {to_cat}', shell=True, capture_output=True)
if p.returncode != 0:
    sys.stderr.write(p.stderr)
    exit()
    
p = subprocess.run(f'rm {to_cat}', shell=True, capture_output=True)
if p.returncode != 0:
    sys.stderr.write(p.stderr)
    exit()


# In[ ]:


sys.stderr.write('completed duplex_subsampler\n')

