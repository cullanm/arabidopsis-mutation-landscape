#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from tqdm import tqdm
from time import time
import argparse
import gtools


# In[ ]:


parser = argparse.ArgumentParser(description='''
Takes a pair of fastq files as input. Adds the first n bp of the read1 and read2 sequence to the end of the first word of the 
sequence id and removes the first m bp from read1 and read2. Designed to be used on libraries made from xGen CS adapters. For 
example, with the default n=3 and m=5:|n
|n
read1 input:|n
@read1_id xxxx|n
CCGATCGAGCTGTAGTCGT|n
|n
read2 input:|n
@read2_id xxxx|n
TTTACGAGGTGGCAGGCGT|n
|n
read1 output:|n
@read1_id_CCGTTT xxxx|n
CGAGCTGTAGTCGT|n
|n
read2 output:|n
@read2_id_CCGTTT xxxx|n
GAGGTGGCAGGCGT
    ''', formatter_class=gtools.EscapeFormatter)
parser.add_argument(dest='fin1', metavar='READ1_FASTQ', type=str, 
                   help='read1 fastq file')
parser.add_argument(dest='fin2', metavar='READ2_FASTQ', type=str, 
                   help='read2 fastq file')
parser.add_argument(dest='fout1', metavar='READ1_OUTPUT', type=str, 
                   help='where to output the modified read1 fastq')
parser.add_argument(dest='fout2', metavar='READ2_OUTPUT', type=str, 
                   help='where to output the modified read2 fastq')
parser.add_argument('-n', '--umi_len', dest='umi_len', metavar='INT', type=int, default=3,
                   help='number of bp to copy from the start of the read to the sequence id (default: 3)')
parser.add_argument('-m', '--adap_len', dest='adap_len', metavar='INT', type=int, default=5,
                   help='number of bp to remove from the start of the read (default: 5)')
parser.add_argument('-e', '--extract_lens', dest='extract_lens', metavar='INT+INT', type=str, default=None,
                   help='single argument for setting -n and -m in format "n+m". Overwrites -n and -m (default: 3+5)')
parser.add_argument('-s', '--separator', dest='separator', metavar='STR', type=str, default='_',
                   help='separator between the initial read id and the UMI (default: _). Note, SAM file format does not allow whitespace in the read id')

try: # run this if in a jupyter notebook
    get_ipython()
    args = parser.parse_args('tests/test_extract_barcode/test_1.fastq tests/test_extract_barcode/test_2.fastq tests/test_extract_barcode/out_1.fastq tests/test_extract_barcode/out_2.fastq'.split()) # used for testing
except: # run this if in a terminal
    args = parser.parse_args()
    
if args.extract_lens is not None:
    args.umi_len = int(args.extract_lens.split('+')[0])
    args.adap_len = int(args.extract_lens.split('+')[-1])


# In[ ]:


fin1 = open(args.fin1, 'r')
fin2 = open(args.fin2, 'r')
fout1 = open(args.fout1, 'w')
fout2 = open(args.fout2, 'w')

j = 0
r1 = ['', '', '', ''] # read1 id, sequence, plus, and quality
r2 = ['', '', '', ''] # read2 id, sequence, plus, and quality
bar = tqdm(position=0, miniters=20000)
# t1 = time()
while True:
#     if j > 1000000:
#         break
    # parse the next 4 lines of each input fastq into r1 and r2
    try:
        for i in range(4):
            r1[i] = next(fin1)
            r2[i] = next(fin2)
    except StopIteration: # if at end of file
        break

#     print(r1)
#     print(r2)
    if r1[0][0] != '@' or r2[0][0] != '@':
        print(f'ERROR: fastq file has unexpected format at line {j * 4}, does not follow a repeating 4 field pattern starting with @')
        break
    r1_words = r1[0].split(' ', 1)
    r2_words = r2[0].split(' ', 1)
    if r1_words[0] != r2_words[0]:
        print('ERROR: read1 and read2 are not in matched order')
        break
    
    r1_bc = r1[1][:args.umi_len]
    r2_bc = r2[1][:args.umi_len]
    r1[0] = r1_words[0] + args.separator + r1_bc + r2_bc + ' ' + r1_words[1]
    r2[0] = r2_words[0] + args.separator + r1_bc + r2_bc + ' ' + r2_words[1]
    r1[1] = r1[1][args.adap_len:]
    r1[3] = r1[3][args.adap_len:]
    r2[1] = r2[1][args.adap_len:]
    r2[3] = r2[3][args.adap_len:]
    
    for i in range(4):
        fout1.write(r1[i])
        fout2.write(r2[i])
    j += 1
    bar.update()
# print(time() - t1)
fin1.close()
fin2.close()
fout1.close()
fout2.close()


# In[ ]:




