#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import pysam
import sys
from tqdm import tqdm
from collections import defaultdict
import os
import subprocess
import multiprocessing
import time
import gtools


# In[ ]:


parser = argparse.ArgumentParser(description='''
Groups read pairs by their UMI and template start and end positions. Adds an integer UMI id tag (UG) and corrected UMI sequence (RX)
tag. The UG tag is unique within the contig. Corrects for sequencing errors by grouping UMIs within n edit distance (hamming distance) and selecting the UMI with the 
greatest count as the corrected UMI. Grouping works by selecting the most common UMI for a given template start and end position, 
grouping all UMIs within n edits from that UMI, then repeating until no UMIs are left. Intended for use with xGen CS adapters. Will discard
any reads not part of a properly aligned pair (flag 2). The output should match that of UMI-tools group but with a few differences: |n
- a different corrected UMI tag (RX)|n
- the ability to group F1R2 and F2R1 reads when UMIs are on each side of the template|n
- marking of both r1 and r2 with the UG and RX tag|n
- can optionally group reads with slightly different template start and end positions
    ''', formatter_class=gtools.EscapeFormatter)
parser.add_argument(dest='fin', metavar='INPUT_BAM', type=str,
                   help='input BAM file, must be coordinate sorted')
parser.add_argument(dest='fout', metavar='OUTPUT_BAM', type=str, 
                   help='output BAM file, will remain coordinate sorted')
parser.add_argument('-s', '--separator', dest='sep', metavar='STR', type=str, default='_',
                   help='separator between the read id and the UMI (default: _)')
parser.add_argument('-f', '--no_flip', dest='flip', action='store_false',
                   help='by default, UMIs of F2R1 reads are flipped from ABCDEF to DEFABC. Setting this flag negates this behavior')
parser.add_argument('-a', '--append', dest='append_umi', metavar='STR', type=str, default='',
                   help='append an arbitrary string to the end of all RX tags, useful when the output will later be \
                   merged with other files which could have UMI conflicts')
parser.add_argument('-e', '--edit_distance', dest='edit_distance', metavar='INT', type=int, default=0,
                   help='max edit (hamming) distance to use when clustering UMIs (default: 0)')
parser.add_argument('-b', '--bp_distance', dest='bp_distance', metavar='INT', type=int, default=0,
                   help='max bp distance the template start and end can differ by when clustering UMIs (default: 0)')
parser.add_argument('-t', '--tmp', dest='tmp', metavar='PREFIX', type=str, default='./',
                   help='prefix to add to temporary files (default: ./)')
parser.add_argument('-@', '--threads', dest='threads', metavar='INT', type=int, default=1,
                   help='number of processors to use, can not utilize more than num contigs + 1 (default=1)')
parser.add_argument('--buffer_size', dest='buffer_size', metavar='INT', type=int, default=20000,
                   help='number of reads kept in memory before flushing to output, adjusting may improve speed (default: 20000)')

try: # run this if in a jupyter notebook
    get_ipython()
    # args = parser.parse_args('tests/test_duplex_caller/duplex_caller_test_case.bam tests/test_duplex_caller/out.tsv -@ 1 --all -t tmp/'.split()) # used for testing
    if os.getcwd()[:8] != '/scratch': # switch to the scratch directory where all the data files are
        os.chdir(f'/scratch/cam02551/{os.getcwd().split("/")[-2]}')
    args = parser.parse_args('-@ 8 -e 0 --append _0 data/align/chandler_KVKCS001D_0_sort.bam data/align/chandler_KVKCS001D_0_umigroup.bam'.split()) # used for testing
except: # run this if in a terminal
    args = parser.parse_args()

if args.fout and '/' in args.fout:
    os.makedirs(os.path.dirname(args.fout), exist_ok=True)

sys.stderr.write(f'Running group_umis.py with arguments {args}\n')


# In[ ]:


# index the input if it hasn't been already
try:
    aln = pysam.AlignmentFile(args.fin, 'rb')
    aln.check_index()
except ValueError:
    sys.stderr.write('No index found for input file, running samtools index\n')
    proc = subprocess.run(f'samtools index -@ {args.threads} {args.fin}', shell=True, capture_output=True)
    if proc.returncode != 0:
        sys.stderr.write('ERROR: failed to run samtools index on input\n')
        sys.stderr.write(proc.stderr)
        sys.exit(1)

# confirm that file is coordinate sorted and has a UMI        
aln = pysam.AlignmentFile(args.fin, 'rb')
last = ('', 0)
for i, read in enumerate(aln.fetch()):
    nex = (read.reference_name, read.reference_start)
    if nex < last:
        sys.stderr.write('ERROR: SAM/BAM file is not coordinate sorted\n')
        exit()
    if len(read.query_name.rsplit(args.sep, 1)) < 2:
        sys.stderr.write(f'ERROR: no separator "{args.sep}" present in read name\n')
    if not set(read.query_name.rsplit(args.sep, 1)[1]).issubset({'A', 'C', 'G', 'T', 'N'}):
        sys.stderr.write(f'WARNING: characters after separator "{args.sep}" do not look like a UMI: {read.query_name.rsplit(args.sep, 1)[1]}\n')
    last = nex
    if i > 1000:
        break
aln.close()


# In[ ]:


# calculates number of edits needed to convert one UMI to another
def hamming_dist(seq1, seq2):
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist


# In[ ]:


def group_umis(fin, fout, chrom, total):
    aln_in = pysam.AlignmentFile(fin, 'rb')
    aln_out = pysam.AlignmentFile(fout, 'wb', template=aln_in)
    read_iter = aln_in.fetch(contig=chrom)

    cur_umi_counts = defaultdict(lambda: defaultdict(lambda: 0)) # cur_umi_counts[fragment end][raw umi] = number of forward reads with the umi
    cur_reads = [] # list of (fragment end, umi, read) tuples
    name_to_umi = dict() # key: read name, value: (umi id, corrected umi)
    cur_id = 0

    pbar = tqdm(miniters=50000, position=0, total=total, desc=chrom, maxinterval=9999)
    read = next(read_iter, None) # get the first read
    pbar.update()
    while not read.is_proper_pair: # skip all reads that aren't proper pairs (flag 2)
            read = next(read_iter, None)
            pbar.update()
    cur_pos = read.reference_start if read is not None else None
    while read is not None:
        # print(read.query_name, read.is_reverse)
        # break
        umi = read.query_name.rsplit(args.sep, 1)[1]
        if args.flip and ((not read.is_reverse and read.is_read2) or (read.is_reverse and read.is_read1)): # flip F2R1 UMIs from ABCDEF to DEFABC
            # print('flipping')
            umi = umi[len(umi) // 2:len(umi)] + umi[:len(umi) // 2]
        tlen = abs(read.template_length)
        cur_reads.append((tlen, umi, read))

        if not read.is_reverse:
            cur_umi_counts[abs(read.template_length)][umi] += 1

        # if read.query_name == 'lh00134:636:22M5YKLT4:7:2231:8378:14340_TGGTTA':
        #     if read.is_reverse:
        #         print('found reverse read')
        #     else:
        #         print('found forward read')
        #     # print(read.is_reverse)
        #     # print(read.reference_start, read.next_reference_start)
        #     print(dict(cur_umi_counts))

        read = next(read_iter, None) # get the next read
        pbar.update()
        while (read is not None) and (not read.is_proper_pair): # skip all reads that aren't proper pairs (flag 2)
            read = next(read_iter, None)
            pbar.update()

        # if read.query_name == 'lh00134:636:22M5YKLT4:7:2231:8378:14340_TGGTTA':
        #     print(f'proper pair? {read.is_proper_pair}')
        #     if read.is_reverse:
        #         print('found reverse read')
        #     else:
        #         print('found forward read')

        # if we advanced to a new start position, flush the reads of the previous one
        if read is None or read.reference_start > cur_pos:
            # print(dict(cur_umi_counts))
            umi_corrector = defaultdict(lambda: dict())
            for tlen in cur_umi_counts:
                sorted_umis = sorted(list(cur_umi_counts[tlen].items()), key=lambda x: x[1], reverse=True)
                # if len(sorted_umis) > 0:
                #     print(f'correcting umis for {cur_pos}+{tlen}: {sorted_umis}')
                for i in range(len(sorted_umis)):
                    if sorted_umis[i][0] not in umi_corrector[tlen]:
                        umi_corrector[tlen][sorted_umis[i][0]] = sorted_umis[i][0]
                        # print(f'umi {sorted_umis[i][0]} is correct')
                    for j in range(i + 1, len(sorted_umis)):
                        if sorted_umis[j][0] in umi_corrector[tlen]:
                            continue
                        if hamming_dist(sorted_umis[i][0], sorted_umis[j][0]) <= args.edit_distance:
                            # print(f'correcting {sorted_umis[j][0]} into {sorted_umis[i][0]}')
                            umi_corrector[tlen][sorted_umis[j][0]] = sorted_umis[i][0]
            # print(f'corrected umis: {dict(umi_corrector)}')
            umi_ids = defaultdict(lambda: dict())
            for tlen, umi, r in cur_reads:
                if r.is_reverse:
                    continue
                # print(tlen, umi)
                cor_umi = umi_corrector[tlen][umi]
                if cor_umi not in umi_ids[tlen]:
                    umi_ids[tlen][cor_umi] = cur_id
                    cur_id += 1
                name_to_umi[r.query_name] = (umi_ids[tlen][cor_umi], cor_umi)

            for tlen, umi, r in cur_reads:
                if not r.is_reverse:
                    cor_umi = umi_corrector[tlen][umi]
                    r.set_tag('UG', umi_ids[tlen][cor_umi])
                    r.set_tag('RX', cor_umi + args.append_umi)
                    # print(f'assigned read {min(r.reference_start, r.next_reference_start)}-{abs(r.tlen)}-{umi} a UG of {umi_ids[tlen][cor_umi]}')
                    aln_out.write(r)
                else:
                    umi_id, cor_umi = name_to_umi[r.query_name]
                    del name_to_umi[r.query_name]
                    r.set_tag('UG', umi_id)
                    r.set_tag('RX', cor_umi + args.append_umi)
                    # print(f'assigned read {min(r.reference_start, r.next_reference_start)}-{abs(r.tlen)}-{umi} a UG of {umi_id}')
                    aln_out.write(r)

            cur_umi_counts = defaultdict(lambda: defaultdict(lambda: 0))
            cur_reads = []
            cur_pos = read.reference_start if read is not None else -1
            # if cur_pos > 2000:
            #     aln_in.close()
            #     aln_out.close()
            #     break


# In[ ]:


# # OLD groups umis for templates which share the same start and end position for one chromosome
# def group_umis(fin, fout, chrom, total): # total is the expected number of reads for drawing the pbar
#     aln_in = pysam.AlignmentFile(fin, 'rb')
#     aln_out = pysam.AlignmentFile(fout, 'wb', template=aln_in)
#     read_buf = [] # list of pysam read objects in order they were read
#     umi_buf = dict() # indexed as umi_buf[(tstart, tlen)][umi] = [read names], dict-dict-list
#     true_umi = dict() # key is read name, value is (corrected_umi, umi_id)
#     cur_id = 1
#     num_read = 0

#     read_iter = aln_in.fetch(contig=chrom)
#     try:
#         read = next(read_iter) # get the first read
#     except StopIteration:
#         sys.stderr.write(f'WARNING: no reads on contig {chrom}\n')
#         aln_in.close()
#         aln_out.close()
#         return
#     pbar = tqdm(miniters=50000, position=0, total=total, desc=chrom)
#     while True:
#     #     print(f'processing read {read.query_name}')
#         # extract the UMI
#         umi = read.query_name.rsplit(args.sep, 1)[1]
#         if args.flip and ((read.is_read1 and read.is_reverse) or (not read.is_read1 and not read.is_reverse)): # flip F2R1 UMIs from ABCDEF to DEFABC
#             umi = umi[len(umi) // 2:len(umi)] + umi[:len(umi) // 2]
#     #         read.query_name = read.query_name[:len(read.query_name) - len(umi)] + umi # FIXME I don't think I need this, but I could add it

#         # add umi and read name to umi_buf
#         coord = (min(read.reference_start, read.next_reference_start), abs(read.template_length))
#         if coord not in umi_buf:
#             umi_buf[coord] = dict()
#         if umi not in umi_buf[coord]:
#             umi_buf[coord][umi] = []
#         umi_buf[coord][umi].append(read.query_name)

#         read_buf.append(read)

#         # get the next properly aligned read, set read to None if at the end of the file
#         while True:
#             read = next(read_iter, None) 
#             if read is None or read.flag & 2:
#                 break

#         # flush reads from the buffer when it gets to large and when the end of the file is reached
#         if len(read_buf) > args.buffer_size or read is None:
#     #         print(read_buf)
#     #         print(umi_buf)
#     #         print(true_umi)
#             for coord in umi_buf:
#                 # skip if we aren't past the template end yet, this is an overestimate. Don't skip if at end of file
#                 if read is not None and coord[0] + coord[1] + 3 > read.reference_start:
#                     continue
#     #             print(umi_buf[frag])
#     #             print(f'processing umis at coords {coord}')
#                 # sort from umi with highest counts to lowest
#                 umis = sorted(list(umi_buf[coord].items()), key=lambda x: len(x[1]), reverse=True) # this is a list of (umi, [read_names]) tuples
#     #             print(umis)
#                 while len(umis) > 0:
#                     mode_umi = umis[0][0]
#     #                 print(f'selected {mode_umi} as mode_umi')
#                     new_umis = []
#                     for i, umi in enumerate(umis):
#                         if hamming_dist(mode_umi, umi[0]) <= args.edit_distance:
#     #                         print(f'merging {umi[0]} with mode_umi')
#                             true_umi.update({x:(mode_umi, cur_id) for x in umi[1]})
#                         else:
#                             new_umis.append(umi)
#                     umis = new_umis
#                     cur_id += 1
#     #         print(true_umi)

#             num_written = 0
#             for r in read_buf:
#                 if r.query_name in true_umi:
#     #                 print(f'writing {r.query_name}')
#                     r.set_tag('UG', true_umi[r.query_name][1]) # add tag for UMI group ID
#                     r.set_tag('RX', true_umi[r.query_name][0] + args.append_umi) # add tag for corrected UMI
#                     aln_out.write(r)
#                     num_written += 1
#                     if r.is_reverse:
#                         del true_umi[r.query_name]
#                 else:
#                     break
#             umi_buf = {x:y for (x, y) in umi_buf.items() if x[0] >= r.reference_start}
#             read_buf = read_buf[num_written:]
#     #         print(f'wrote {num_written} reads, {len(read_buf)} remain. umi_buf size of {len(umi_buf)}, true_umi size of {len(true_umi)}')
#             pbar.update(num_written)
#     #         if num_written == 0:
#     #             sys.stderr.out('ERROR: exceeded read buffer size, this could be a bug or --buffer_size may need to be larger')
#     #             sys.exit(1)

#         # break if at the end of the file
#         if read is None:
#             del umi_buf
#             break

#     aln_in.close()
#     aln_out.close()
#     sys.stderr.write(f'Done grouping {chrom}\n')


# In[ ]:


# group UMIs which vary by a few bp in start or end position, OLD and probably really slow
def group_neighbors(fin, fout, chrom, total):
    aln_in = pysam.AlignmentFile(fin, 'rb')
    aln_out = pysam.AlignmentFile(fout, 'wb', template=aln_in)
    read_buf = [] # list of pysam read objects in order they were read
    cur_pos = -1
    group_buf = dict() # indexed as umi_buf[(tstart, tlen)] = {(umi, umi_id)} dict-set
    new_id = dict() # used to convert umis being grouped to the new one key:umi_id value:umi_id to group the key into
    id_to_umi = dict()
    # merged_umis = dict() # inverse of new_umi, finds the UMIs being grouped into a umi key:umi_id value:(true_umi, umi_id being grouped into the key)
    group_counts = defaultdict(lambda: 0) # key:umi_id value:number of reads with that id 

    read_iter = aln_in.fetch(until_eof=True)
    try:
        read = next(read_iter) # get the first read
    except StopIteration:
        sys.stderr.write(f'WARNING: no reads on contig {chrom}\n')
        aln_in.close()
        aln_out.close()
        return
    pbar = tqdm(miniters=20000, position=0, total=total, desc=chrom)
    while True:
    #     print(f'processing read {read.query_name}')
        # extract the UMI

        # add umi and id to buffer
        coord = (min(read.reference_start, read.next_reference_start), abs(read.template_length))
        if coord not in group_buf:
            group_buf[coord] = set()
        group_buf[(min(read.reference_start, read.next_reference_start), abs(read.template_length))].add((read.get_tag('RX'), read.get_tag('UG')))
        group_counts[read.get_tag('UG')] += 1
        id_to_umi[read.get_tag('UG')] = read.get_tag('RX')

        read_buf.append(read)

        read = next(read_iter, None) # get the next read, set read to None if at the end of the file

        # flush reads from the buffer when it gets to large and when the end of the file is reached
        if len(read_buf) > args.buffer_size or read is None:
    #         print(read_buf)
    #         print(umi_buf)
    #         print(true_umi)
            for coord in group_buf: # for each (templat start, template len) in the buffer
    #             print(f'joining neighbors around {coord}')
                if read is not None and coord[0] + coord[1] + 3 + args.bp_distance > read.reference_start: # skip if we could still find more to group
                    continue
                for central_umi in group_buf[coord]: # for each (umi, umi_id) at the coord position
                    for i in range(coord[0] - args.bp_distance, coord[0] + args.bp_distance + 1): # vary template start
                        for j in range(coord[1] - args.bp_distance, coord[1] + args.bp_distance + 1): # vary template len
                            if i == coord[0] and j == coord[1]: # ignore identical template coord, these should already be grouped
                                continue
                            if (i, j) not in group_buf: # ignore template positions with no reads
                                continue
                            for neighbor_umi in group_buf[(i, j)]: # for each (umi, umi_id) at the varied template position
                                if hamming_dist(central_umi[0], neighbor_umi[0]) <= args.edit_distance: # if the umis are similar, group
    #                                 print(f'grouping {central_umi} {coord} and {neighbor_umi} {(i, j)}')
                                    # group the umi with fewer reads into the one with more
                                    if group_counts[central_umi[1]] > group_counts[neighbor_umi[1]]:
                                        new_id[neighbor_umi[1]] = central_umi[1]
                                    elif group_counts[central_umi[1]] < group_counts[neighbor_umi[1]]:
                                        new_id[central_umi[1]] = neighbor_umi[1]
                                    else:
                                        new_id[min(central_umi[1], neighbor_umi[1])] = max(central_umi[1], neighbor_umi[1])


            num_written = 0
            for r in read_buf:
                # stop writing once we find a read that might still need grouping, don't stop if at end of file
                if read is not None and r.reference_start + r.tlen + 3 + args.bp_distance > read.reference_start:
                    break
    #             print(f'writing {r.query_name}')
    #             print(new_id)
                new = r.get_tag('UG')
                while new in new_id:
                    new = new_id[new]
                r.set_tag('UG', new) # add tag for UMI group ID
                r.set_tag('RX', id_to_umi[new]) # add tag for corrected UMI
                aln_out.write(r)
                num_written += 1
            to_remove = [y for (x, y) in group_buf.items() if x[0] + x[1] + args.bp_distance < r.reference_start]
            for s in to_remove:
                for u in s:
                    del id_to_umi[u[1]]
                    new_id.pop(u[1], None)
            group_buf = {x:y for (x, y) in group_buf.items() if x[0] + x[1] + args.bp_distance >= r.reference_start}
            read_buf = read_buf[num_written:]
            pbar.update(num_written)
#             print(f'wrote {num_written} reads, {len(read_buf)} remain. group_buf:{len(group_buf)}, id_to_umi:{len(id_to_umi)}, new_id:{len(new_id)}')


        # break if at the end of the file
        if read is None:
            break

    aln_in.close()
    aln_out.close()
    sys.stderr.write(f'Done grouping {chrom}\n')


# In[ ]:


if __name__ == '__main__':
    aln = pysam.AlignmentFile(args.fin, 'rb')
    chroms = {x.contig:x.mapped for x in aln.get_index_statistics()}
    chroms = {x:chroms[x] for x in sorted(chroms.keys())}

    # if the tmp file directory doesn't exist, make it
    if not os.path.exists(args.tmp.rsplit('/', 1)[0]):
        os.makedirs(args.tmp.rsplit('/', 1)[0])

    # get file names for temporary files
    tmp_files = {chrom:args.tmp + args.fout.replace('/', '_') + '_' + chrom + '.bam' for chrom in chroms}
    sys.stderr.write(f'using temporary files {", ".join(tmp_files.values())}\n')

    # start worker processes from pool
    sys.stderr.write('Grouping umis with the same template start and end\n')
    if args.threads > 1:
        with multiprocessing.Pool(processes=args.threads - 1) as pool:
            processes = []
            for chrom in chroms: # for each chromosome, add a process
                processes.append(pool.apply_async(group_umis, (args.fin, tmp_files[chrom], chrom, chroms[chrom])))
            pool.close()
            pool.join()
        sys.stderr.flush()
        time.sleep(2)
        sys.stderr.write('Checking for errors during runtime\n')
        # check for errors
        for p in processes:
            p.get()
    else:
        for chrom in chroms:
            group_umis(args.fin, tmp_files[chrom], chrom, chroms[chrom])


# In[ ]:


if __name__ == '__main__' and args.bp_distance > 0:
    tmp_nei_files = {chrom:args.tmp + args.fout.rsplit('.', 1)[0] + '_' + chrom + '_nei.bam' for chrom in chroms}
    sys.stderr.write(f'using another set of temporary files {", ".join(tmp_nei_files.values())}\n')

    # start worker processes from pool
    sys.stderr.write(f'Grouping umis off by up to {args.bp_distance}bp in template start and/or end\n')
    if args.threads > 1:
        with multiprocessing.Pool(processes=args.threads - 1) as pool:
            processes = []
            for chrom in chroms: # for each chromosome, add a process
                processes.append(pool.apply_async(group_neighbors, (tmp_files[chrom], tmp_nei_files[chrom], chrom, chroms[chrom])))
            pool.close()
            pool.join()
        sys.stderr.flush()
        time.sleep(2)
        sys.stderr.write('Checking for errors during runtime\n')
        # check for errors
        for p in processes:
            p.get()
    else:
        for chrom in chroms:
            group_neighbors(tmp_files[chrom], tmp_nei_files[chrom], chrom, chroms[chrom])

    sys.stderr.write(f'removing initial temporary files\n')
    for file in tmp_files.values():
        try:
            os.remove(file)
        except FileNotFoundError:
            pass

    tmp_files = tmp_nei_files


# In[ ]:


sys.stderr.write('Concatenating temporary files into final output\n')
proc = subprocess.run(f'samtools cat -@ {args.threads} -o {args.fout} {" ".join(tmp_files.values())}', shell=True, capture_output=True)
if proc.returncode != 0:
    sys.stderr.write('ERROR: failed to run samtools cat on temporary files\n')
    sys.stderr.write(proc.stderr)
    sys.exit(1)

# sys.stderr.write('Removing any remaining temporary files\n')
# for file in tmp_files.values():
#     try:
#         os.remove(file)
#     except FileNotFoundError:
#         pass

sys.stderr.write(f'completed group_umis.py on {args.fin}\n')


# In[ ]:




