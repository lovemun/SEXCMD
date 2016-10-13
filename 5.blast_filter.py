#!/usr/bin/env python

import sys, os, string, operator, glob, gzip, math, time

args = sys.argv
if len(args) != 4:
    print 'Usage:', args[0], '.inp.fasta .blastm9 .out.fasta'
    sys.exit()

fastafilename = args[1]
blastfilename = args[2]
outfafilename = args[3]

sexdict = {'chrX':0, 'chrY':0, 'chrZ':0, 'chrW':0}

fadict = {}
fahead = ''
for lines in open(fastafilename, 'r'):
    if lines.startswith('>'):
        fahead = lines[1:].strip()
        fadict[fahead] = ''
    else:
        fadict[fahead] = lines.strip()

blastdict = {}
for lines in open(blastfilename, 'r'):
    if lines.startswith('#'): continue
    linelist = lines.rstrip('\n').split('\t')
    Query_id, Subject_id, pIdentity, alignment_length, mismatches, gap_openings, q_start, q_end, \
        s_start, s_end, e_value, bit_score = linelist
    pCoverage = float(int(q_end)-int(q_start)-1)/float(len(fadict[Query_id]))
    # 70% Identity and 70% Coverage for Matching Criteria
    if pCoverage >= 0.7 and float(pIdentity) >= 70.:
        if int(alignment_length) >= 50: # min match 50bp
            blastdict.setdefault(Query_id, {}).setdefault(Subject_id, []).append(linelist)

blastlist = blastdict.keys()
blastlist.sort()

finaldict = {}
for Query_id in blastlist:
    chrList = blastdict[Query_id].keys()
    chrList.sort()
    if len(chrList) != 2: continue # only two sex chromosomes
    if len([ix for ix in chrList if sexdict.has_key(ix)]) != 2: continue # only two sex chromosomes
    iskip = 0
    for sexchr in sexdict.keys():
        if not blastdict[Query_id].has_key(sexchr): continue
        if len(blastdict[Query_id][sexchr]) > 1: iskip = 1
    if iskip:
        #print >> sys.stderr, 'more than two hits in each sex chromosome', blastdict[Query_id]
        continue
    finaldict[Query_id] = 0

# identity filter
fMinPercentIdentity = 0.05
fMaxPercentIdentity = 0.9

finalList = [ix for ix in finaldict.keys() if ix.startswith('chrX') or ix.startswith('chrZ')]
finalList.sort()

outfile = open(outfafilename, 'w')

for XID in finalList:
    # chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20
    i1, i2, i3, i4 = XID.split('|')
    YID = i2 + '|' + i1 + '|' + i3 + '|' + i4
    outfile.write('>%s\n' % XID)
    outfile.write('%s\n' % fadict[XID])
    outfile.write('>%s\n' % YID)
    outfile.write('%s\n' % fadict[YID])

outfile.close()


"""
# BLASTN 2.2.26 [Sep-21-2011]
# Query: chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20
# Database: hg19
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. e
nd, e-value, bit score
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	100.00	197	0	0	1	197	6995261	
6995457	1e-106	 391
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	100.00	18	0	0	179	196	49395733
	49395750	7.1	36.2
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	92.31	26	2	0	3	28	94038939
	94038964	7.1	36.2
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	95.45	22	1	0	106	127	12000778
3	120007762	7.1	36.2
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	95.45	22	1	0	106	127	12006417
6	120064155	7.1	36.2
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	95.45	22	1	0	106	127	12006903
7	120069016	7.1	36.2
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	95.45	22	1	0	106	127	12007389
8	120073877	7.1	36.2
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	95.45	22	1	0	106	127	12007875
8	120078737	7.1	36.2
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	95.45	22	1	0	106	127	12008361
9	120083598	7.1	36.2
chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20	chrX	95.45	22	1	0	106	127	12008847
9	120088458	7.1	36.2

[deepreds@my hg19]$ head sex_marker*fasta
>chrX:6995261-6995457;+|chrY:17460549-17460745;+|197|20
CTTCTCCATAGCAGGAGGGGGAGAGAACCTCTTGGCACAAGCTAGGAAGATGTCCGGGTCTGGCTTGCCATGCTGCACTTCGGGGTCATCTCCCAGCACAATGTGGGAAAACAAGCTGAAGAACTCCTTGTGGCGGCTTGTCTTCATATCGAACGACGCGGACCCCGAGCTGGTGGCCAGTGCAAAGGGGATGCCAT
"""

