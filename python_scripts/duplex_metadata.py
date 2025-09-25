#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import os
import gtools
import pysam
import matplotlib.pyplot as plt
from tqdm import tqdm
from collections import defaultdict
import random
import numpy as np
import statistics as stat
import matplotlib
import pandas as pd
import matplotlib.colors as colors
import multiprocessing
import pickle
plt.rcParams['svg.fonttype'] = 'none'
import itertools
import argparse


# In[ ]:


parser = argparse.ArgumentParser(description='''
Outputs metadata statistics and figures on a list of NanoSeq BAM files. Outputs {prefix}metadata.tsv with columns:|n
sample rep sample_read_pairs_whole_genome sample_conflict_rate rep_read_pairs rep_frags rep_frags_both_strands rep_frags_call_overlap rep_frags_call_whole
- sample: name of BAM file|n
- rep: name of replicate extracted from the UMI tag (optional)|n
- sample_read_pairs_whole_genome: total number of concordant read pairs in the sample. All other statistics are for only the region specified with --region|n
- sample_conflict_rate: the estimated fraction of fragments in the sample which will not be callable due to another fragment at the same alignment position 
with sufficient reads|n
- rep_read_pairs: concordant read pairs in the replicate (in --region)|n
- rep_frags: number of fragments in the replicate (in --region)|n
- rep_frags_both_strands: number of fragments in the replicate (in --region) with both strands sequenced |n
- rep_frags_call_overlap: number of fragments in the replicate (in --region) with enough coverage to call in the read1-read2 overlap (>=B/2 read pairs per 
strand and >=C/2 read pairs total)|n
- rep_frags_call_whole: number of fragments in the replicate (in --region) with enough coverage to call the whole fragment (>=B read pairs per strand and 
>=C read pairs total)|n
|n
Outputs the following figures:|n
- {prefix}conflict_rates.svg: scatterplot showing the relationship between fragment count and conflict rate as well as the line used to impute conflict 
rate for samples with one replicate|n
- {prefix}subsampling_efficiency.svg: used to estimate dilution efficiency. Replicates libraries have their reads subsampled to various levels and sequencing
efficiency (callable fragments / aligned read pair) is calculated for each level. An optimally diluted library will plateau at 1, while an oversequenced library
will peak at a subsample level < 1. Note: callable fragments here does not take into account mapping quality or blacklists.|n
- {prefix}reads_per_molecule.svg: histogram of the number of reads per molecule.|n
|n
Lastly, a dictionary of the conflict rates is written to stdout for copying into the duplex_mutation_figures notebook.
    ''', formatter_class=gtools.EscapeFormatter)
parser.add_argument('-i', '--input', dest='input', metavar='BAM', type=str, nargs='+', required=True,
                   help='input BAM files of duplex-seq libraries')
parser.add_argument('-o', '--output', dest='output', metavar='PREFIX', type=str, required=True,
                   help='prefix to prepend the output metadata table and figures')
parser.add_argument('-u', '--umi', dest='umi', metavar='TAG', type=str,
                   help='BAM tag containing the UMI (default: no UMI)')
parser.add_argument('-r', '--replicate', dest='replicate', metavar='SEPARATOR,INDEX', type=str,
                   help='should be set if a single BAM file contains multiple replicate samples. These replicates will be \
                   treated separately in the metadata output and are necessary for calculating an estimated conflict \
                   rate for a BAM file. The replicate name is extracted from the UMI tag by splitting \
                   on SEPARATOR and taking the substring at INDEX. e.g. with "--umi RX --sample_extract _,1", \
                   an RX tag of "CACGTC_library1" will extract "library1" as the replicate name. Pass an \
                   empty SEPARATOR to indicate that the UMI is the replicate name (default: treat each \
                   BAM as one replicate)')
parser.add_argument('-g', '--region', dest='region', metavar='CHROM:START-END', type=str,
                   help='consider only reads with start positions within this region. Significantly reduces runtime \
                   (default: use whole genome)')
parser.add_argument('-s', '--subsample_levels', dest='sub_levels', metavar='FLOAT', type=float, nargs='+', default=[0.2, 0.4, 0.6, 0.8, 1],
                   help='subsample the library to these fractions of reads when visualizing dilution efficiency \
                   (default: 0.2 0.4 0.6 0.8 1)')
parser.add_argument('-B', '--min_strand_cov', dest='min_strand_cov', metavar='INT', type=int, default=2,
                   help='filter cutoff for calling variants, used for determining fragments with enough coverage to call: \
                   minimum reads per strand overlapping variant (default: 2)')
parser.add_argument('-C', '--min_total_cov', dest='min_total_cov', metavar='INT', type=int, default=6,
                   help='filter cutoff for calling variants, used for determining fragments with enough coverage to call: \
                   minimum total reads overalpping variant (default: 6)')
parser.add_argument('-N', '--min_frac_final', dest='min_frac_final', metavar='FLOAT', type=float, default=0.76,
                   help='filter cutoff for calling variants, used in calculating conflict rate: final minimum fraction \
                   of supporting reads per strand (default: 0.76)')
parser.add_argument('-@', '--threads', dest='threads', metavar='INT', type=int, default=1,
                   help='multiprocessing threads for parallel processing of BAM files (default: 1)')
parser.add_argument('--tmp', dest='tmp', metavar='PREFIX', type=str, default='tmp/',
                   help='prefix for temporary files (default: ./)')

try: # run this if in a jupyter notebook
    get_ipython()
    print('Determined code is running in Jupyter')
    if os.getcwd()[:8] != '/scratch': # switch to the scratch directory where all the data files are
        os.chdir(f'/scratch/cam02551/{os.getcwd().split("/")[-2]}')
    args = parser.parse_args('-i data/align/big_Col-0-1_merged.bam data/align/big_Col-0-2_merged.bam -o tmp/metadata_test_ -g Chr1:0-100000 -u RG -r ,0 --tmp tmp/'.split()) # used for testing
except: # run this if in a terminal
    args = parser.parse_args()

sys.stderr.write('Running duplex_metadata.py with arguments:\n' + '\n'.join([f'{key}={val}' for key, val in vars(args).items()]) + '\n')

if args.output and '/' in args.output:
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
if '/' in args.tmp:
    os.makedirs(os.path.dirname(args.tmp), exist_ok=True)


# In[ ]:


os.path.dirname('asdf')


# In[ ]:


if args.region:
    contig = args.region.split(':')[0]
    pos_start = int(args.region.split(':')[1].split('-')[0])
    pos_end = int(args.region.split('-')[1])

if args.replicate:
    if not args.umi:
        sys.stderr.write('ERROR: must specificy --umi if specifying --replicate')
        sys.exit()
    rep_sep = args.replicate.split(',')[0]
    rep_index = int(args.replicate.split(',')[1])


# # Extract fragment information from the BAM files

# In[ ]:


# creates a pickled dictionary with key=replicate value=[(chrom, start, tlen, umi, + strand pairs, - strand pairs)]
def get_frags(file, out, pbar_i, sub=2):
    aln = pysam.AlignmentFile(file, 'rb')
    total_reads = sum([x.total for x in aln.get_index_statistics()])
    genome_size = sum([aln.get_reference_length(x.contig) for x in aln.get_index_statistics()])
    pbar_total = int((pos_end - pos_start) / genome_size * total_reads) if args.region else total_reads
    fetch_iter = aln.fetch(contig=contig, start=pos_start, end=pos_end) if args.region else aln.fetch()

    fragments = defaultdict(lambda: [0, 0])
    rep_reads = defaultdict(lambda: 0)
    for read in tqdm(fetch_iter, position=pbar_i, total=pbar_total, miniters=200000, maxinterval=999):
        frag_start = read.reference_start
        tlen = abs(read.template_length) # this is non-inclusive
        if args.umi:
            umi = read.get_tag(args.umi)
        else:
            umi = ''
        if args.replicate:
            rep = umi.split(rep_sep)[rep_index] if len(rep_sep) > 0 else umi
        else:
            rep = ''

        rep_reads[rep] += 1

        # thow out read pairs that don't align as expected, ignore read 2 so we don't double count pairs
        if read.is_unmapped or read.mate_is_unmapped or not (read.flag & 2) or read.is_reverse or read.mapping_quality < 10:
            continue

        if random.random() > sub:
            continue

        frag_strand = '+' if read.is_read1 else '-'

        if frag_strand == '+':
            fragments[(read.reference_name, frag_start, tlen, umi, rep)][0] += 1
        else:
            fragments[(read.reference_name, frag_start, tlen, umi, rep)][1] += 1

    fragments = dict(fragments)
    frags = {rep:[] for rep in sorted(rep_reads)}
    for frag in fragments:
        frags[frag[4]].append((frag[0], frag[1], frag[2], frag[3], fragments[frag][0], fragments[frag][1]))

    with open(out, 'wb') as fout:
        pickle.dump(frags, fout)

    aln.close()

    # tqdm.write(f'{os.path.basename(file)}: {sum(rep_reads.values())} reads{" in " + args.region if args.region else ""}', file=sys.stdout)
    # for rep in rep_reads:
    #     tqdm.write(f'\treplicate {rep}: {rep_reads[rep]} reads{" in " + args.region if args.region else ""}', file=sys.stdout)


# In[ ]:


if __name__ == '__main__':
    tqdm.write('counting reads and fragments in BAM files\n', file=sys.stderr)
    pool = multiprocessing.Pool(processes=8)
    processes = []
    for i, f in enumerate(args.input):
        processes.append(pool.apply_async(get_frags, (f, f'{args.tmp}duplex_metadata_{os.path.basename(f)}_frags.pkl', i)))
    pool.close()
    pool.join()
    # check for errors
    for p in processes:
        p.get()


# In[ ]:


tqdm.write('loading temporary fragment files\n', file=sys.stderr)
frags = dict() # key: sample; value: [(chrom, frag_start, tlen, umi, + reads, - reads)]
for file_name in tqdm(args.input, position=0):
    with open(f'{args.tmp}duplex_metadata_{os.path.basename(file_name)}_frags.pkl', 'rb') as f:
        l = pickle.load(f)
        frags[os.path.basename(file_name)] = l
    os.remove(f'{args.tmp}duplex_metadata_{os.path.basename(file_name)}_frags.pkl')


# # Plot reads per fragment

# In[ ]:


# plot total coverage of each fragment to estimate the number of viable fragments in the library
fig, axs = plt.subplots(sum([len(frags[s]) for s in frags]), figsize=(4, 0.5 * sum([len(frags[s]) for s in frags])), sharex=True)
i = 0
for sample in tqdm(frags):
    for rep in frags[sample]:
        axs[i].hist([x[4] + x[5] for x in frags[sample][rep]], bins=range(30))
        axs[i].hist([x[4] + x[5] for x in frags[sample][rep] if x[4] > 0 and x[5] > 0], bins=range(30))
        axs[i].set_title(f'{sample} {rep}', x=0.98, y=0.4, ha='right', size=8)
        i += 1
axs[-1].set_xlabel('# read pairs covering fragment')
axs[len(axs) // 2].set_ylabel('# molecules')
axs[0].legend(['Reads in 1 strand', 'Reads in both strands'], loc='lower center', bbox_to_anchor=(0.5, 1))
fig.suptitle('Read coverage per molecule', y=0.94)
fig.savefig(f'{args.output}reads_per_molecule.svg', dpi=300, bbox_inches='tight')


# # Calculate conflict rate

# In[ ]:


tqdm.write('calculating conflict rates\n', file=sys.stderr)
conflict_rates = dict()
for i, sample in tqdm(enumerate(frags), total=len(frags), position=0):

    # for each replicate, randomly sample half as many fragments from other replicates to compute a "half" conflict rate 
    # (the true conflict rate is calculated as 1 - (1-half conflict rate)^2)
    # the reason I calculate a "half" rate is that some samples only have 2 replicates, and if one has more fragments than the other, its 
    # "full" conflict rate can't be calculated, so I use the half as a close approximation
    rates = []
    rep_frag_counts = []
    if len(frags[sample]) == 1:
        conflict_rates[sample] = 0
        continue

    for rep in frags[sample]:
        rep_frags = frags[sample][rep] # fragments from this rep
        other_frags = list(itertools.chain.from_iterable([frags[sample][r] for r in frags[sample] if r != rep])) # fragments from other reps

        # randomly select half as many fragments from other_frags as are in rep_frags
        ratio = len(rep_frags) / len(other_frags)
        if ratio > 2:
            continue
        other_subsample = {(x[0], x[1], x[2]):(x[4], x[5]) for x in other_frags if random.random() < 0.5 * ratio}

        # count the number of conflicts which make a fragment uncallable
        uncallable_conflicts = 0
        for frag in rep_frags:
            if (frag[0], frag[1], frag[2]) in other_subsample: # if there's a conflict
                # if the fragment doesn't contribute at least 76% of the top and bottom strand pairs, the conflict will be uncallable
                if frag[4] / max(1, (frag[4] + other_subsample[(frag[0], frag[1], frag[2])][0])) < args.min_frac_final or \
                   frag[5] / max(1, (frag[5] + other_subsample[(frag[0], frag[1], frag[2])][1])) < args.min_frac_final:
                    uncallable_conflicts += 1
        conflict_rate = 1 - (1 - (uncallable_conflicts / len(rep_frags))) ** 2
        rates.append(conflict_rate)

        rep_frag_counts.append(sum([x[4] > 0 and x[5] > 0 for x in rep_frags])) # estimated number of callable fragments

    # calculate the sample conflict rate as the average across replicates, weighted by the number of fragments in the replicate
    conflict_rates[sample] = sum([rates[i] * rep_frag_counts[i] / sum(rep_frag_counts) for i in range(len(rates))])


# In[ ]:


# make a best fit correlation between fragment number and conflict rate (only use samples where conflict rate could be estimated)
counts = []
rates = []
for sample in frags:
    if conflict_rates[sample] == 0:
        continue
    counts.append(sum([len(frags[sample][r]) for r in frags[sample]]))
    rates.append(conflict_rates[sample])
counts += [0] * 10000 # my tiny brain can't figure out how to fix the intercept at 0, so I'm just adding a bunch of points at (0, 0)
rates += [0] * 10000
slope, intercept = np.polyfit(counts, rates, 1) if len(counts) > 10000 else (0, 0)
line_x = np.linspace(0, 5000000, 1000)
line_y = [x * slope + intercept for x in line_x]


# In[ ]:


imputed_rates = conflict_rates.copy()
for sample in imputed_rates:
    if imputed_rates[sample] == 0:
        imputed_rates[sample] = sum([len(frags[sample][r]) for r in frags[sample]]) * slope + intercept


# In[ ]:


fig, ax = plt.subplots()
ax.scatter([sum([len(frags[s][r]) for r in frags[s]]) for s in imputed_rates], imputed_rates.values(), color=['tab:blue' if conflict_rates[s] == 0 else 'tab:orange' for s in imputed_rates])
xlim, ylim = ax.get_xlim(), ax.get_ylim()
ax.plot(line_x, line_y)
ax.set_xlim(xlim), ax.set_ylim(ylim)
ax.set_xlabel(f'Fragment number in {args.region}' if args.region else 'Fragment number')
ax.set_ylabel('Conflict rate')
ax.set_title('Conflict rate and fragment number')
s_blue = plt.scatter([], [], c='tab:blue')
s_orange = plt.scatter([], [], c='tab:orange')
ax.legend([s_orange, s_blue], ['calculated', 'imputed'])
fig.savefig(f'{args.output}conflict_rate.svg', dpi=300, bbox_inches='tight')


# # Create metadata tsv

# In[ ]:


tqdm.write('writing metadata tsv\n', file=sys.stderr)
metadata = []
for f in args.input:
    sample = os.path.basename(f)
    aln = pysam.AlignmentFile(f, 'rb')
    sample_pairs = sum([x.total for x in aln.get_index_statistics()]) // 2
    aln.close()
    for rep in frags[sample]:
        rep_pairs = sum([x[4] + x[5] for x in frags[sample][rep]])
        rep_frags = len(frags[sample][rep])
        rep_both_strands = sum([x[4] > 0 and x[5] > 0 for x in frags[sample][rep]])
        rep_call_overlap = sum([x[4] >= args.min_strand_cov // 2 and x[5] >= args.min_strand_cov // 2 and x[4] + x[5] >= args.min_total_cov // 2 for x in frags[sample][rep]])
        rep_call_whole = sum([x[4] >= args.min_strand_cov and x[5] >= args.min_strand_cov and x[4] + x[5] >= args.min_total_cov for x in frags[sample][rep]])
        metadata.append((sample, rep, sample_pairs, imputed_rates[sample], rep_pairs, rep_frags, rep_both_strands, rep_call_overlap, rep_call_whole))
df_meta = pd.DataFrame(metadata, columns='sample rep sample_read_pairs_whole_genome sample_conflict_rate rep_read_pairs rep_frags rep_frags_both_strands rep_frags_call_overlap rep_frags_call_whole'.split())
df_meta.to_csv(f'{args.output}metadata.tsv', sep='\t', index=False)


# # Subsample to determine sequencing efficiency

# In[ ]:


sub_levels = sorted(args.sub_levels)


# In[ ]:


if __name__ == '__main__':
    tqdm.write('counting reads and fragments while subsampling BAM files\n', file=sys.stderr)
    pool = multiprocessing.Pool(processes=8)
    processes = []
    for i, f in enumerate(args.input):
        for j, sub_level in enumerate(sub_levels):
            processes.append(pool.apply_async(get_frags, (f, f'{args.tmp}duplex_metadata_{os.path.basename(f)}_{sub_level}_frags.pkl', i * len(args.input) + j, sub_level)))
    pool.close()
    pool.join()

    # check for errors
    for p in processes:
        p.get()


# In[ ]:


tqdm.write('loading temporary subsampled fragment files\n', file=sys.stderr)
sub_frags = defaultdict(lambda: dict())
for i, file_name in tqdm(enumerate(args.input)):
    for sub_level in sub_levels:
        with open(f'{args.tmp}duplex_metadata_{os.path.basename(file_name)}_{sub_level}_frags.pkl', 'rb') as f:
            frags = pickle.load(f)
            for replicate in sorted(list(frags.keys())):
                sub_frags[os.path.basename(file_name) + ' ' + replicate][sub_level] = frags[replicate]
        os.remove(f'{args.tmp}duplex_metadata_{os.path.basename(file_name)}_{sub_level}_frags.pkl')


# In[ ]:


subset_frags = {x:sub_frags[x] for x in list(sub_frags)}

cmap = plt.cm.tab20
fig, axs = plt.subplots(3, figsize=(len(subset_frags), 6), sharex=True, sharey=True)
for i, req_strand, req_total in zip(range(3), [1, 1, 2], [2, 3, 4]):
    seq_eff_overlap = []
    seq_eff_all = []
    for sample in subset_frags:
        for sub_level in sub_levels:
            total_reads = sum([x[4] + x[5] for x in subset_frags[sample][sub_level]])
            callable_frags_overlap = len([x for x in subset_frags[sample][sub_level] if x[4] >= req_strand and x[5] >= req_strand and x[4] + x[5] >= req_total])
            callable_frags_all = len([x for x in subset_frags[sample][sub_level] if x[4] >= 2 * req_strand and x[5] >= 2 * req_strand and x[4] + x[5] >= 2 * req_total])
            seq_eff_overlap.append(callable_frags_overlap / total_reads)
            seq_eff_all.append(callable_frags_all / total_reads)
    axs[i].bar(range(len(seq_eff_overlap)), seq_eff_overlap, color=[cmap((x // len(sub_levels)) % 2) for x in range(len(seq_eff_overlap))])
    axs[i].bar(range(len(seq_eff_all)), seq_eff_all, color=[cmap((x // len(sub_levels)) % 2 + 2) for x in range(len(seq_eff_all))])
    axs[i].set_title(f'Need >={req_strand * 2} reads per strand and >={req_total * 2} total to call', y=0.7)
axs[1].set_ylabel('Callable fragmens per concordant read pair')
axs[-1].set_xticks([len(sub_levels) // 2 + x * len(sub_levels) for x in range(len(subset_frags))])
axs[-1].set_xticklabels(subset_frags.keys(), rotation=45, ha='right')
axs[0].legend(['only read1/2 overlap is callable', 'all of fragment is callable'], bbox_to_anchor=(1, 1), loc='upper left')
blank = plt.Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
axs[1].legend([blank] * 2, ['subsample levels:', ', '.join(map(str, sub_levels))] + sub_levels, bbox_to_anchor=(1, 1), loc='upper left')
fig.suptitle('Dilution efficiency of subsampled libraries', y=0.95)
fig.savefig(f'{args.output}subsampling_efficiency.svg', dpi=300, bbox_inches='tight')


# In[ ]:


tqdm.write('conflict rates:\n', file=sys.stdout)
tqdm.write(str(imputed_rates), file=sys.stdout)
tqdm.write(f'completed duplex_metadata.py on {args.input}\n')

