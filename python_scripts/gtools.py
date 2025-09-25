#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import subprocess
from tqdm import tqdm
from itertools import chain
from multiprocessing import Process
import numpy as np
import pysam
# import matplotlib.pyplot as plt


# In[ ]:


def load_genome(file):
    genome = dict()
    with open(file, 'r') as f:
        ls = f.read()
        ls = ls.split('>')[1:]
        for l in ls:
            genome[l.split(None, 1)[0]] = l.split('\n', 1)[1].replace('\n', '')
    return genome


# In[ ]:


def running_notebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


# In[ ]:


# Load a VCF file as a pandas dataframe with format matching that of duplex_caller. 
# POS: changed to 0-based
# REF & ALT: converts indels to asterisk format, and unpacks records containing multiple variants at the same position to separate rows
# e.g. the variant Chr1 10 TA C,TCAA,T,TAG,TAAA,TAAAA will be converted to:
#                  Chr1 9  T  C
#                  Chr1 10 A     C
#                  Chr1 11 *      AA
#                  Chr1 10 A         *
#                  Chr1 11 *             G
#                  Chr1 11 *                 AA
#                  Chr1 11 *                      AAA
# 
# INFO: unpacked to separate columns for each info item in the header
# FORMAT & SAMPLEs: each sample column is unpacked into a separate column for each format field in the header
def load_vcf(file, descriptions=False, nrecords=None):
    fin = open(file, 'r')
    l = next(fin, None)

    infos = dict()
    formats = dict()
    while l is not None:
        if l[0] == '#':
            if l[1] != '#':
                samples = l[1:].strip().split()[9:]
            elif l[:6] == '##INFO':
                infos[l.split('ID=')[1].split(',')[0]] = (l.split('Number=')[1].split(','), l.split('Description="')[1].split('"')[0])
            elif l[:8] == '##FORMAT':
                formats[l.split('ID=')[1].split(',')[0]] = (l.split('Number=')[1].split(','), l.split('Description="')[1].split('"')[0])
        else:
            break
        l = next(fin, None)

    header = 'chrom pos id ref alt qual fil'.split() + list(infos.keys())
    for sample in samples:
        header += [f'{sample}_{form}' for form in formats]

    variants = []
    records_read = 0
    while l is not None:
        var_base = {h:None for h in header}
        raw = l.strip().split('\t')

        # add the static fields
        var_base['chrom'] = raw[0]
        var_base['pos'] = raw[1]
        var_base['id'] = raw[2]
        var_base['ref'] = raw[3]
        var_base['alt'] = raw[4]
        var_base['qual'] = raw[5]
        var_base['fil'] = raw[6]

        # add the INFO fields
        for x in raw[7].split(';'):
            x = x.split('=')
            if len(x) == 1:
                val = True
            elif infos[x[0]] == 'A':
                val = x[1].split(',')
            else:
                val = x[1]
            var_base[x[0]] = val

        # add the FORMAT fields
        forms = raw[8].split(':')
        for sample_idx, sample_field in enumerate(raw[9:]):
            for format_idx, value in enumerate(sample_field.split(':')):
                var_base[f'{samples[sample_idx]}_{forms[format_idx]}'] = value

        start_pos = int(var_base['pos'])
        ref = var_base['ref']
        alts = var_base['alt'].split(',')
        vars_to_add = []
        for alt_idx, alt in enumerate(alts):
            for i in range(min(len(ref), len(alt))): # iterate through the positions until we run out of bases in ref or alt
                if ref[i] != alt[i]: # output all mismatches as snps
                    v = var_base.copy()
                    v['pos'] = start_pos - 1 + i
                    v['ref'] = ref[i]
                    v['alt'] = alt[i]
                    for info in [x for x in infos if infos[x][0] == 'A']:
                        v[info] = var_base[info][alt_idx]
                    vars_to_add.append(v)
            if len(ref) == len(alt): # if there's no indel here
                pass
            elif len(ref) > len(alt): # if a deletion, the starting bases will be lined up, so the extra bases at the end of ref are what's deleted
                v = var_base.copy()
                v['pos'] = start_pos - 1 + min(len(ref), len(alt))
                v['ref'] = ref[len(alt):]
                v['alt'] = '*'
                for info in [x for x in infos if infos[x][0] == 'A']:
                    v[info] = var_base[info][alt_idx]
                vars_to_add.append(v)
            else: # in an insertion, the extra bases at the end of alt are what's inserted
                v = var_base.copy()
                v['pos'] = start_pos - 1 + min(len(ref), len(alt))
                v['ref'] = '*'
                v['alt'] = alt[len(ref):]
                for info in [x for x in infos if infos[x][0] == 'A']:
                    v[info] = var_base[info][alt_idx]
                vars_to_add.append(v)
        vars_to_add.sort(key=lambda v: v['pos'])
        variants += [list(v.values()) for v in vars_to_add]

        l = next(fin, None)
        records_read += 1
        if nrecords is not None and records_read == nrecords:
            break

    fin.close()

    df = pd.DataFrame(variants, columns=header)
    df.pos = pd.to_numeric(df.pos)
    if descriptions:
        return df, infos, formats
    else:
        return df


# In[ ]:


# import sys
# import os
# sys.path.append(os.getcwd() + '/../python_scripts') # this lets us import files in python_scripts (like gtools)
# import gtools
# if os.getcwd()[:8] != '/scratch': # switch to the scratch directory where all the data files are
#     os.chdir(f'/scratch/cam02551/{os.getcwd().split("/")[-2]}')

# df = load_vcf('data/variant/mutect/big_Col-0-1_filtered.vcf')


# In[ ]:


def load_chrom_sizes(index):
    chrom_sizes = dict()
    with open(index, 'r') as f:
        for l in f:
            split = l[:-1].split('\t')
            chrom_sizes[split[0]] = int(split[1])
    return chrom_sizes


# In[ ]:


def load_coverage_bedgraph(bg, index):
    chrom_sizes = load_chrom_sizes(index)
    coverage = {chrom:np.zeros(chrom_sizes[chrom], dtype='i') for chrom in chrom_sizes}
    with open(bg, 'r') as f:
        for l in f:
            split = l[:-1].split('\t')
            coverage[split[0]][int(split[1]):int(split[2])] = int(float(split[3]))
    return coverage


# In[ ]:


def run_commands(commands):
    for command in tqdm(commands):
        # print('Running: ' + command)
        comp = subprocess.run(command, shell=True, capture_output=True)
        if comp.returncode != 0:
            print(comp)
    return comp.returncode


# In[ ]:


def add_duplex_color_tag(fin, fout):
    if fin[-4:] == '.sam':
        aln_in = pysam.AlignmentFile(fin, 'r')
        aln_out = pysam.AlignmentFile(fout, 'w', template=aln_in)
    elif fin[-4:] == '.bam':
        aln_in = pysam.AlignmentFile(fin, 'rb')
        aln_out = pysam.AlignmentFile(fout, 'wb', template=aln_in)
    else:
        print('ERROR: unknown file type')

    for read in aln_in.fetch():
        frag_start = min(read.reference_start, read.next_reference_start)
        frag_len = abs(read.template_length)
        if read.has_tag('RX'):
            frag_umi = read.get_tag('RX')
        elif read.has_tag('RG'):
            frag_umi = read.get_tag('RG')
        else:
            frag_umi = ''
        read.set_tag('co', f'{frag_start}_{frag_len}_{frag_umi}')
        aln_out.write(read)


# In[ ]:


# Uses SAMtools to pick out alignment file reads around positions listed in a dataframe
# df must contain a column named chom and pos
# files is a list of sorted, indexed BAM file names
# window_size: will output reads within window_size bp of each variant
# run: if false, will print commands rather than running them
# max_reads: max number of reads to 
def sam_view_df(df, files, outdir='data/align/views/', window_size=100, run=True, max_reads=5000, color=True):
    positions = list(set(zip(df.chrom, df.pos))) # store positions as (chrom, pos) and drop duplicates
    positions.sort(key=lambda x: x[1]) # sort by position
    positions.sort(key=lambda x: x[0]) # sort by chromosome, this results in proper genome position order

    # group positions that are withing 2 * window_size from each other
    last_pos = ('asdfasdf', -100000)
    groups = [] 
    for pos in positions:
        if pos[0] == last_pos[0] and pos[1] - last_pos[1] <= 2 * window_size: # if positions are too close
            groups[-1].append(pos)
        else:
            groups.append([pos])
        last_pos = pos

    # convert groups of positions into windows for viewing
    windows = [] # (chrom, start, end)
    for group in groups:
        windows.append((group[0][0], max(0, group[0][1] - window_size), group[-1][1] + window_size))

    # generate and run/print samtools view and sort commands
    folder = files[0].rsplit('/', 1)[0]
    commands_list = []
    tmp_files = [f'sam_view_df_tmp_{window[0]}_{window[1]}_{window[2]}.sam' for window in windows]
    for file in files:
        commands = []
        for i, window in enumerate(windows):
            commands.append(f'samtools view {file} {window[0]}:{window[1]}-{window[2]} | head -n {max_reads} > {tmp_files[i]}')
        commands.append(f'samtools view -H {file} > sam_view_df_tmp_header.sam')
        commands.append(f'cat sam_view_df_tmp_header.sam {" ".join(tmp_files)} | samtools sort -o {outdir}view_{file.rsplit("/", 1)[-1].rsplit(".", 1)[0]}.sam')
        commands_list.append(commands)
    if run:
        for command_list in commands_list:
            for command in command_list:
                # print('Running: ' + command)
                comp = subprocess.run(command, shell=True, capture_output=True)
                if comp.returncode != 0:
                    print(comp)
        comp = subprocess.run(f'rm {" ".join(tmp_files)} sam_view_df_tmp_header.sam', shell=True, capture_output=True)
        if comp.returncode != 0:
            print(comp)

        if color:
            for file in files:
                add_duplex_color_tag(f'{outdir}view_{file.rsplit("/", 1)[-1].rsplit(".", 1)[0]}.sam', 'tmp.sam')
                subprocess.run(f'mv tmp.sam {outdir}view_{file.rsplit("/", 1)[-1].rsplit(".", 1)[0]}.sam', shell=True, capture_output=True)
    else:
        for commands in commands_list:
            print(' &\n'.join(commands[:-1]))
            print('wait')
            print(commands[-1])
        print(f'rm {" ".join(tmp_files)} sam_view_df_tmp_header.sam')


# In[ ]:


# iterable for streaming multiple bam files in position sorted order, functionally the same as reading a file generated by samtools merge
class BamStreamer:
    aln_files = []
    streams = []
    cur_reads = []
    read_starts = []
    def __init__(self, bam_files, chrom):
        self.aln_files = [pysam.AlignmentFile(f, 'rb') for f in bam_files]
        self.streams = [self.aln_files[i].fetch(contig=chrom) for i in range(len(bam_files))]
        self.cur_reads = [next(s) for s in self.streams]
        self.read_starts = [r.reference_start for r in self.cur_reads]
    def __iter__(self):
        return self
    def __next__(self):
        min_idx = min(range(len(self.read_starts)), key=self.read_starts.__getitem__)
        to_return = self.cur_reads[min_idx]
        self.cur_reads[min_idx] = next(self.streams[min_idx])
        self.read_starts[min_idx] = self.cur_reads[min_idx].reference_start
        return min_idx, to_return
    def close(self):
        for f in self.aln_files:
            f.close()


# In[ ]:


# rainbow order: green, 
# the colors at index 9 and 10 
__base_colors = ['#08843E', '#6ABD45', '#CEC82B', '#FDB714', '#F6891F', '#D55127', '#CE2127', '#C61E6F', '#A63493', '#6D3A96', '#4C479D', '#3960AD', '#1B99D5', '#29AD95', '#8F5423']
__dark_colors = ['#09542B', '#508D40', '#9D9936', '#CB942B', '#C56F29', '#A54023', '#9F2428', '#961D55', '#762468', '#4E2667', '#34326D', '#28467D', '#005C99', '#1F7D6C', '#5F3819']
__light_colors = ['#7CC57C', '#9ECC42', '#FAEE3B', '#FDCE57', '#F5A04E', '#F06941', '#F04D60', '#E05CA2', '#B86DAD', '#8866AC', '#696AB1', '#5D80C1', '#76B3E2', '#7ECCBF', '#BA795E']

# the colors at index 2&3 (yellow-green and goldenrod), 5&6 (dark orange and red) and 9&10 (dark purple and dark blue) may be indistinguishable for certain types of colorblindness, 
# thus, excluding the last 3 colors makes the colormap more colorblind friendly
# the red and brown may also be difficult to distinguish for red-blind protanopia, but it depended on the simulation

__reorder = [0, 1, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14]

base_cmap = [__base_colors[i] for i in __reorder]
dark_cmap = [__dark_colors[i] for i in __reorder]
light_cmap = [__light_colors[i] for i in __reorder]

full_cmap = []
for i in __reorder:
    full_cmap.append(__dark_colors[i])
    full_cmap.append(__base_colors[i])
    full_cmap.append(__light_colors[i])


# In[ ]:


# import matplotlib
# %matplotlib inline
# fig, ax = matplotlib.pyplot.subplots()
# for i, color in enumerate(base_cmap):
#     ax.add_patch(matplotlib.patches.Rectangle((0 + i, 0), 1, 1, color=color))
# for i, color in enumerate(dark_cmap):
#     ax.add_patch(matplotlib.patches.Rectangle((0 + i, 1), 1, 1, color=color))
# for i, color in enumerate(light_cmap):
#     ax.add_patch(matplotlib.patches.Rectangle((0 + i, -1), 1, 1, color=color))
# ax.set_xlim(0, len(base_cmap))
# ax.set_ylim(-1, 2)


# In[ ]:


# formatter to pass argparse formatter_class which 
import argparse
class EscapeFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        import textwrap
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = textwrap.fill(paragraph.replace('|t', '\t'), width, initial_indent=indent, subsequent_indent=indent) + '\n'
            multiline_text += formatted_paragraph
        return multiline_text 
    # def _fill_text(self, text, width, indent):
    #     text = self._whitespace_matcher.sub(' ', text).strip()
    #     import textwrap
    #     return textwrap.fill(text, width,
    #                          initial_indent=indent,
    #                          subsequent_indent=indent).replace('|n', '\n').replace('|t', '\t')


# In[ ]:




