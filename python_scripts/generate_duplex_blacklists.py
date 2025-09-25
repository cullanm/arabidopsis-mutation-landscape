#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import os
from tqdm import tqdm
import random
from collections import defaultdict
import re
import sys
sys.path.append(os.getcwd() + '/../python_scripts')
import gtools
import argparse
import subprocess
import pysam
import multiprocessing
import math


# In[ ]:


parser = argparse.ArgumentParser(description='''
Creates blacklist files for filtering variants of duplex-seq libraries. Takes the reference genome and a set of BAM files for coverage calculations as input. 
Requires GenMap (https://github.com/cpockrandt/genmap) to be installed. Produces four .npy arrays for each chromosome: |n
    - {prefix}{chrom}_coverage_percentile.npy: the coverage percentile of each site in the genome, e.g. 0.01 means 1% of the genome has a lower or equal coverage|n
    - {prefix}{chrom}_duplicate_distance.npy: the minimum hamming distance between a kmer at the genomic position of length --fragment_length and any other position in the genome.
        e.g. 2 means there's another kmer in the genome where only 2 bases differ. If no kmer is discovered, a value of 10 is assigned.|n
    - {prefix}{chrom}_poly_at_length.npy: length of the longest poly-A/T repeat starting at the position or 1bp away. Only considers lengths >=4. e.g.|n
|t        sequence:  ACTGAAAAAAAAGTCCCA|n
|t        blacklist: 000880000008800000|n
    - {prefix}{chrom}_poly_repeat_length.npy: length of the longest poly-di/trinucleotide repeat starting at the position of 1bp away. Only considers length >=3. e.g.|n
|t        sequence:  CCCGAGAGAGAGACCTC|n
|t        blacklist: 00554000000455000
    ''', formatter_class=gtools.EscapeFormatter)
parser.add_argument('-f', '--fasta', required=True, dest='fasta', metavar='FASTA', type=str,
                   help='genome fasta')
parser.add_argument('-b', '--bedgraphs', required=True, dest='bedgraphs', metavar='BAMs', type=str, nargs='+',
                   help='one or more coverage bedgraph files used to calculate coverage percentile per genomic site. These can \
                   be generated with "bedtools genomecov -bg -ibam {BAM}"')
parser.add_argument('-o', '--output', dest='prefix', metavar='PREFIX', type=str, default='./',
                   help='prefix to prepend to blacklist outputs (default: ./)')
parser.add_argument('-l', '--fragment_length', dest='fragment_length', metavar='INT', type=int, default=150,
                   help='length of DNA fragments. Used for calculating the duplicate distance blacklist and passed to genmap -k argument (default: 150)')
parser.add_argument('-m', '--max_hamming', dest='max_hamming', metavar='INT', type=int, default=3,
                    help='maximum hamming distance to check with genmap for making the duplicate distance blacklist. Increasing this increases runtime. Setting \
                    to a negative value will skip running genmap and the blacklist will be all 10s (default: 3)')
parser.add_argument('-@', '--threads', dest='threads', metavar='INT', type=int, default=1,
                    help='number of threads to use when calculating genome coverage (default: 1)')
parser.add_argument('-t', '--tmp', dest='tmp_prefix', metavar='PREFIX', type=str, default='tmp/',
                   help='prefix to add to temporary files. Will make any directories which do not exist (default: tmp/)')
# parser.add_argument('-@', '--threads', dest='threads', metavar='INT', type=int, default=1,
#                    help='number of threads. Each chromosome can be processed by a separate thread. (default: 1)')

try: # run this if in a jupyter notebook
    get_ipython()
    print('Determined code is running in Jupyter')
    if os.getcwd()[:8] != '/scratch': # switch to the scratch directory where all the data files are
        os.chdir(f'/scratch/cam02551/{os.getcwd().split("/")[-2]}')
    args = parser.parse_args(' -@ 8 --fasta data/ref/ref.fa --bams data/align/big_Col-0-1_merged.bam data/align/big_Col-0-1_merged.bam data/align/big_Col-0-1_merged.bam data/align/big_Col-0-1_merged.bam data/align/big_Col-0-1_merged.bam data/align/big_Col-0-1_merged.bam data/align/big_Col-0-7_merged.bam data/align/big_Col-0-8_merged.bam --output data/blacklist/new_arabidopsis/'.split()) # used for testing
except: # run this if in a terminal
    args = parser.parse_args()

sys.stderr.write('Running generate_duplex_blacklists.py with arguments:\n' + '\n'.join([f'{key}={val}' for key, val in vars(args).items()]) + '\n')

no_match_dist = 10
if args.max_hamming >= 10:
    no_match_dist = args.max_hamming + 1
    print(f'WARNING: --max_hamming is >=10, changing the value reported for no kmer detected to {args.max_hamming + 1}')

if args.prefix and '/' in args.prefix:
    os.makedirs(os.path.dirname(args.prefix), exist_ok=True)
if args.tmp_prefix and '/' in args.tmp_prefix:
    os.makedirs(os.path.dirname(args.tmp_prefix), exist_ok=True)


# In[ ]:


def run_shell(command):
    proc = subprocess.run(command, shell=True, capture_output=True)
    if proc.returncode != 0:
        sys.stderr.write(f'ERROR: failed to run "{command}"\n')
        sys.stderr.write(str(proc.stderr))
        sys.exit(1)


# In[ ]:


genome = gtools.load_genome(args.fasta)
chrom_sizes = {chrom:len(genome[chrom]) for chrom in genome}


# # Duplicate distance blacklists

# In[ ]:


if args.max_hamming >= 0: # only run this if necessary
    sys.stderr.write('calculating duplicate distance blacklist\n')

    sys.stderr.write('making a genmap index of the reference fasta\n')
    index_dir = f'{args.tmp_prefix}generate_duplex_blacklists_genmap_index_{args.fasta.replace("/", "_")}'
    run_shell(f'genmap index -F {args.fasta} -I {index_dir}')

    # for each hamming distance (allowable mismatches), run genmap map
    duplicate_dist = {chrom:np.full(chrom_sizes[chrom], no_match_dist, dtype=np.byte) for chrom in chrom_sizes}
    for i in range(args.max_hamming, -1, -1): 
        sys.stderr.write(f'searching for kmers with max of {i} errors\n')
        bedgraph_prefix = f'{args.tmp_prefix}generate_duplex_blacklists_{args.fasta.replace("/", "_")}_{args.fragment_length}_{i}'
        run_shell(f'genmap map -K {args.fragment_length} -bg -E {i} -I {index_dir} -O {bedgraph_prefix}')

        # open the genmap bedgraph and add its info to duplicate_dist
        with open(f'{bedgraph_prefix}.bedgraph', 'r') as f:
            for l in f:
                fields = l.split()
                # set values within the bedgraph interval. If occurences is 1 (unique) keep value the same, else set value to i
                duplicate_dist[fields[0]][int(fields[1]):int(fields[2])] = duplicate_dist[fields[0]][int(fields[1]):int(fields[2])] * (fields[3] == '1') + i * (fields[3]  != '1')
        run_shell(f'rm {bedgraph_prefix}.bedgraph')
    run_shell(f'rm -r {index_dir}')

    for i in range(args.max_hamming + 1):
        sys.stderr.write(f'{sum([np.count_nonzero(duplicate_dist[chrom] == i) for chrom in chrom_sizes]) / sum(chrom_sizes.values()) * 100}% of the genome has a (120, {i}) duplicate\n')

    sys.stderr.write(f'writing duplicate distance blacklists to {args.prefix}{{chrom}}_duplicate_distance.npy\n')
    for chrom in duplicate_dist:
        np.save(f'{args.prefix}{chrom}_duplicate_distance.npy', duplicate_dist[chrom])
    del duplicate_dist
else:
    sys.stderr.write(f'--max_hamming is <0, creating duplicate distance blacklists with a {no_match_dist} for every site\n')
    for chrom in chrom_sizes:
        np.save(f'{args.prefix}{chrom}_duplicate_distance.npy', np.full(chrom_sizes[chrom], no_match_dist, dtype=np.byte))


# # Coverage blacklists

# In[ ]:


# def genomecov(fin, fout):
#     run_shell(f'bedtools genomecov -bg -ibam {fin} > {fout}')


# In[ ]:


# convert BAM files to bedgraphs
# if __name__ == '__main__':
#     # start worker processes from pool
#     sys.stderr.write('converting BAM files to coverage bedgraphs using bedtools genomecov\n')
#     tmp_files = [f'{args.tmp_prefix}generate_duplex_blacklists_coverage_{args.fasta.replace("/", "_")}_{i}.bg' for i in range(len(args.bedgraphs))]
#     with multiprocessing.Pool(processes=max(1, args.threads)) as pool:
#         processes = []
#         for i, f in enumerate(args.bedgraphs): # for each chromosome, add a process
#             processes.append(pool.apply_async(genomecov, (f, tmp_files[i])))
#         pool.close()
#         pool.join()
#     sys.stderr.flush()
#     time.sleep(2)
#     # check for errors
#     for p in processes:
#         p.get()


# In[ ]:


sys.stderr.write('loading coverage bedgraphs\n')
coverage = {chrom:np.zeros(chrom_sizes[chrom], dtype=int) for chrom in chrom_sizes}
for bg in args.bedgraphs:
    with open(bg, 'r') as f:
        for l in f:
            fields = l.split()
            if 'e' in fields[3]: # if number is in scientific notation, python needs to convert to int using a float intermediate
                c = int(float(fields[3]))
            else:
                c = int(fields[3])
            coverage[fields[0]][int(fields[1]):int(fields[2])] += c
    # run_shell(f'rm {tmp_file}')


# In[ ]:


# randomly sample 10M positions of the genome and get their coverage
sys.stderr.write('calculating conversion of coverage to percentile\n')
sample_size = 1000000
sampled_covs = []
for chrom in chrom_sizes:
    k = int(sample_size * (chrom_sizes[chrom] / sum(chrom_sizes.values())))
    for i in range(k):
        sampled_covs.append(coverage[chrom][random.randrange(chrom_sizes[chrom])])

# count the number of positions with each coverage value
maxcov = int(np.percentile(sampled_covs, 99.9))
cov_counts = np.zeros(maxcov + 1)
for cov in sampled_covs:
    cov_counts[min(maxcov, cov)] += 1

# # count the number of positions with coverage less than n for n=(0, max coverage)
# cov_less_than = np.zeros(maxcov)
# for i in range(len(cov_less_than)):
#     cov_less_than[i] = sum(cov_less_thha)

# make an array that takes a coverage value and outputs the percentile
percentile = np.zeros(maxcov)
for i in range(len(percentile)):
    percentile[i] = sum(cov_counts[:i]) / len(sampled_covs)


# In[ ]:


sys.stderr.write(f'coverage\tpercentile\n')
next_print = 0
for i in range(len(percentile)):
    if percentile[i] > next_print:
        sys.stderr.write(f'{i}\t\t{percentile[i]}\n')
        next_print = math.ceil(percentile[i])


# In[ ]:


# convert each coverage value in the genome into a percentile
sys.stderr.write('converting coverage values to percentiles\n')
cov_percentile = dict() # cov_percentile[chrom][pos] = percent of genome with lesser coverage
for chrom in tqdm(chrom_sizes):
    # apply function to coverage array to convert coverage to percentile
    cov_percentile[chrom] = np.vectorize(lambda x: percentile[x] if x < len(percentile) else percentile[-1])(coverage[chrom])
    del coverage[chrom]


# In[ ]:


sys.stderr.write(f'writing coverage percentile blacklists to {args.prefix}{{chrom}}_coverage_percentile.npy\n')
for chrom in cov_percentile:
    np.save(f'{args.prefix}{chrom}_coverage_percentile.npy', cov_percentile[chrom])
del cov_percentile


# # Poly A/T and poly di/tri-nucleotide repeat blacklists

# In[ ]:


sys.stderr.write('calculating poly A/T and poly di/tri-nucleotide repeat blacklists\n')


poly_at_len = {chrom:np.zeros(chrom_sizes[chrom], dtype='B') for chrom in chrom_sizes} # unsigned byte max is 256
    # poly_at_len[chrom][pos] = length of longest poly-A or poly-T which starts or ends at the position (or the position is 1bp outside of)
poly_rep_len = {chrom:np.zeros(chrom_sizes[chrom], dtype='B') for chrom in chrom_sizes}
    # poly_rep_len[chrom][pos] = length (in number of repeats) of longest di/tri-nuc repeat which starts or ends at the position (or the position is 1bp outside of)

# get regex iterables for each possible nucleotide repeat
match_iters = defaultdict(lambda: [])
for chrom in chrom_sizes:
    match_iters[chrom].append((1, re.finditer('AAAA+', genome[chrom]))) # poly A, the "1" here is the length of repeat, will be "2" for di and "3" for tri
    match_iters[chrom].append((1, re.finditer('TTTT+', genome[chrom]))) # poly T
    for x in 'ACGT':
        for y in 'ACGT':
            if x != y: # two nucleotides must be different to be a dinucleotide repeat
                match_iters[chrom].append((2, re.finditer('(' + x + y + '){3,}', genome[chrom]))) # dinucleotide repeats
            for z in 'ACGT':
                if x != y or x != z or y != z: # at least one pair of nucleotides must be different to be a trinucleotide repeat
                    match_iters[chrom].append((3, re.finditer('(' + x + y + z + '){3,}', genome[chrom]))) # trinucleotide repeats

# iterate through each regex iterable and add the length of the repeat to the array
for chrom in match_iters:
    for replen, it in tqdm(match_iters[chrom], desc=chrom): # for each regex search (e.g. CT repeats)
        for match in it:
            start = match.span()[0]
            end = match.span()[1]
            poly_len = min(256, (end - start) // replen) # unsigned bytes have a max value of 256
            if replen == 1:
                poly_at_len[chrom][start - 1] = max(poly_at_len[chrom][start - 1], poly_len)
                poly_at_len[chrom][start] = max(poly_at_len[chrom][start], poly_len)
                poly_at_len[chrom][end - 1] = max(poly_at_len[chrom][end - 1], poly_len)
                if end != chrom_sizes[chrom]: # this happens if the chrom ends with a poly repeat
                    poly_at_len[chrom][end] = max(poly_at_len[chrom][end], poly_len)
            else:
                poly_rep_len[chrom][start - 1] = max(poly_rep_len[chrom][start - 1], poly_len)
                poly_rep_len[chrom][start] = max(poly_rep_len[chrom][start], poly_len)
                poly_rep_len[chrom][end - 1] = max(poly_rep_len[chrom][end - 1], poly_len)
                if end != chrom_sizes[chrom]:
                    poly_rep_len[chrom][end] = max(poly_rep_len[chrom][end], poly_len)

for i in range(4, 10):
    sys.stderr.write(f'blacklisted {sum(np.count_nonzero(poly_at_len[chrom] == i) for chrom in chrom_sizes) / sum(chrom_sizes.values()) * 100}% for length {i} A/T nucleotide repeats\n')
for i in range(3, 10):
    sys.stderr.write(f'blacklisted {sum(np.count_nonzero(poly_rep_len[chrom] == i) for chrom in chrom_sizes) / sum(chrom_sizes.values()) * 100}% for length {i} di/tri nucleotide repeats\n')


# In[ ]:


sys.stderr.write(f'writing poly-A/T blacklists to {args.prefix}{{chrom}}_poly_at_length.npy\n')
sys.stderr.write(f'writing poly-di/trinucleotide repeat blacklists to {args.prefix}{{chrom}}_poly_repeat_length.npy\n')
for chrom in chrom_sizes:
    np.save(f'{args.prefix}{chrom}_poly_at_length.npy', poly_at_len[chrom])
    np.save(f'{args.prefix}{chrom}_poly_repeat_length.npy', poly_rep_len[chrom])


# In[ ]:


sys.stderr.write(f'completed generate_duplex_blacklists.py on fasta {args.fasta}\n')

