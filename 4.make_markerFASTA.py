#!/usr/bin/env python

#chrX	6995261	6995457	+	CTTCTCCATAGCAGGAGGGGGAGAGAACCTCTTGGCACAAGCTAGGAAGATGTCCGGGTCTGGCTTGCCATGCTGCACTTCGGGGTCATCTCCCAGCACAATGTGGGAAAACAAGCTGAAGAACTCCTTGTGGCGGCTTGTCTTCATATCGAACGACGCGGACCCCGAGCTGGTGGCCAGTGCAAAGGGGATGCCAT	6995269,6995270,6995278,6995297,6995315,6995331,6995341,6995361,6995363,6995383,6995408,6995411,6995417,6995420,6995424,6995425,6995427,6995429,6995441,6995456	ENST00000412827,ENST00000381077,ENST00000424830,ENST00000540122,ENST00000486446
#chrY	17460549	17460745	+	CTTCTCCACGGCAGGAGTGGGAGAGAACCTCTTGGCGCAAGCTAGGAAGATGTCGGGGTCTGGCTTGCCACGCTGCACTTTGGGGTCATCTCCCAGCACAGTCTGGGAAAACAAGCTGAAGATCTCCTTGTGGCGGCTTGTCTTCATCTCAAACGATGCAGACTTCAACCTGGTGGCCAGCGCAAAGGGGATGCCTT	17460557,17460558,17460566,17460585,17460603,17460619,17460629,17460649,17460651,17460671,17460696,17460699,17460705,17460708,17460712,17460713,17460715,17460717,17460729,17460744	ENST00000421058
#chrX	7023764	7023879	+	CTGTGCCGCCTCTAATGCCTTCTTACCCATAACCAGGGACTTTACATCCCAGCTGTATTTCTTGTCATAGCGATTACATATTTCTTGAAACACCACTGAATACAGCCGTTCAGTAT	7023770,7023782,7023783,7023808,7023826,7023827,7023834,7023838,7023850,7023863,7023871	ENST00000486446,ENST00000424830,ENST00000498474,ENST00000412827,ENST00000540122,ENST00000381077
#chrY	17493607	17493722	+	CTGTGCTGCCTCTAATGCAATCTTACCCATAACCAGGGACTTTATATCCCAGCTGTATTTCTCATCATAGTGATCACATATTTCTTCAAACACCACTGAGTACAGCCATTCAGTAT	17493613,17493625,17493626,17493651,17493669,17493670,17493677,17493681,17493693,17493706,17493714	ENST00000421058

import sys,os,string

args = sys.argv
if len(args) != 2:
    print 'Usage : '+ args[0] + ' sex_marker_filtered.txt'
    sys.exit()

infile = open(args[1], 'r')
outfile = open(args[1][:-4] + '.fasta', 'w')
while 1:
    first_line = infile.readline()
    if not first_line: break
    first_line_component = first_line.strip().split('\t')
    second_line = infile.readline()
    second_line_component = second_line.strip().split('\t')
    
    if len(first_line_component[4]) != len(second_line_component[4]): continue
    else:
        cnt = 0
        for idx in range(len(first_line_component[4])):
            if first_line_component[4][idx] != second_line_component[4][idx]:
                cnt = cnt + 1
        outfile.write('>%s:%s-%s;%s|%s:%s-%s;%s|%d|%d\n'%(first_line_component[0], first_line_component[1], first_line_component[2], first_line_component[3], second_line_component[0], second_line_component[1], second_line_component[2], second_line_component[3], len(first_line_component[4]), cnt))
        outfile.write(first_line_component[4]+'\n')
        outfile.write('>%s:%s-%s;%s|%s:%s-%s;%s|%d|%d\n'%(second_line_component[0], second_line_component[1], second_line_component[2], second_line_component[3], first_line_component[0], first_line_component[1], first_line_component[2], first_line_component[3], len(second_line_component[4]), cnt))
        outfile.write(second_line_component[4]+'\n')

outfile.close()
