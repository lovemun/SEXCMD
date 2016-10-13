#!/usr/bin/env python

import sys, os
from itertools import groupby

args = sys.argv
if len(args) != 3:
    print 'Usage :' + args[0] + ' input.fasta chr'
    sys.exit()

infile = open(args[1], 'r')
CHR = 'chr'+str(args[2])
outfile = open(args[1][:-5]+CHR+'.fasta', 'w')

chr_dict = {}
seq = ''
for is_header, group in groupby(infile, lambda x:x[0] == '>'):
    if is_header:
        chr = group.next()[1:].strip()
    else:
        glist = map(lambda x:x.strip(), group)
        seq = ''.join(glist)
    chr_dict[chr] = seq

outfile.write('>' + CHR + '\n')
for idx in xrange(0, len(chr_dict[CHR]), 50):
    outfile.write(chr_dict[CHR][idx:idx+50] + '\n')
outfile.close()

#>chr1
#CATGACacttttgaacaatttcatggggtggtttggacggggaaaggcta
#gccccataagtgaagttttgccagggcacatcccgggattcctaggatcc
#ccatatcatgggattgcttgtattattgaaaaatggggtcccctacggga
#tgcaggaaatgttcaggaatccccagcccgcctggctgagtatttcagta
#gtctaaacaaaaacgtactcactacacgagaacgacagttagctgcctcc
#acggtggcgtggccgctgctgatggccttgaaccggctaaatataactga
#aggggttctcactgatgaaaatcaactattacgtgaccgtgtggaagaac
#tggaaagacaggttgcaattttgagagggaaaaaacctaaccctataatt
#ccggtacaaaaaatccgaaacatctctctagaaataactggggatcctat
