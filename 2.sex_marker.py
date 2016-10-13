#!/usr/bin/env python

import sys
import os.path
import subprocess

class openMAF:
	def __init__(self, file):
		self.file = file
		self.reader = open(file, 'r')
	def __iter__(self):
		return self
	def readParagraph(self):
		paragraph = ''
		for line in self.reader:
			if line.rstrip('\n') == '': break
			if line[0] == '#': continue
			paragraph += line
		return paragraph
	def next(self):
		paragraph = self.readParagraph()
		if paragraph == '':
			raise StopIteration
		return paragraph
	def close(self):
		self.reader.close()

class Alignment:
	# 1-based inclusive
	def __init__(self, chromosome, start, end, strand, sequence):
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.strand = strand
		self.sequence = sequence
		self.index2position = {}
		self.position2index = {}
		if self.strand == '+':
			indexList = range(len(self.sequence))
		if self.strand == '-':
			indexList = range(len(self.sequence))[::-1]	# reverse
		position = self.start
		for index in indexList:
			if not self.sequence[index] == '-':		# not match(indel)
				self.index2position[index] = position
				self.position2index[position] = index
				position += 1
		self.index2genes = {}
	# Annotate genes using tabix
	def annotateGene(self, gtfFile):
		pipe = subprocess.Popen("tabix %s.gz %s:%d-%d" % (gtfFile, self.chromosome, self.start, self.end), shell=True, stdout=subprocess.PIPE).stdout
		for line in pipe:
			tokens = line.rstrip('\n').split('\t')
			if tokens[2] == "exon":
				for attribute in tokens[8].strip().rstrip(';').split(';'):
					attribute = attribute.strip()
					(key, value) = attribute.split(' ', 1)
					value = value.strip()
					if value[0] == '"' and value[-1] == '"':
						value = value[1:-1]
					if key == "gene_id":
						(start, end) = map(int, [tokens[3], tokens[4]])
						(startIndex, endIndex) = sorted([self.position2index[max(start, self.start)], self.position2index[min(end, self.end)]])
						endIndex += 1
						for index in range(startIndex, endIndex):
							self.index2genes.setdefault(index, set()).add(value)
		pipe.close()

def printMarkers(alignments, startIndex, endIndex, positionsList, genesListList):
	positionsList = [[x[i] for x in positionsList] for i in range(len(alignments))]
	genesListList = [[x[i] for x in genesListList] for i in range(len(alignments))]
	for alignment, positions, genesList in zip(alignments, positionsList, genesListList):
		genes = set([gene for genes in genesList for gene in genes])
		sequence = alignment.sequence[startIndex:endIndex]
		# start, end: 1-based inclusive
		(start, end) = sorted([alignment.index2position[startIndex], alignment.index2position[endIndex - 1]])
		print '\t'.join([alignment.chromosome, str(start), str(end), alignment.strand, sequence, ','.join(map(str, positions)), ','.join(genes)])

def getStartIndex(alignments, startIndex, minimumStartIndex):
	while startIndex > minimumStartIndex:
		if any([alignment.sequence[startIndex - 1] == '-' for alignment in alignments]):
			break
		if any([not alignment.index2genes.has_key(startIndex - 1) for alignment in alignments]):
			break
		startIndex -= 1
	return startIndex

def getEndIndex(alignments, endIndex, maximumEndIndex):
	while endIndex < maximumEndIndex:
		if any([alignment.sequence[endIndex] == '-' for alignment in alignments]):
			break
		if any([not alignment.index2genes.has_key(endIndex) for alignment in alignments]):
			break
		endIndex += 1
	return endIndex

# Input
(inputFile, gtfFile, maximumDistance) = sys.argv[1:]
if not os.path.isfile("%s.gz.tbi" % gtfFile):
	subprocess.call("grep -v '^#' %s | sort --field-separator='\t' -k1,1 -k4,4n -k5,5n | bgzip > %s.gz" % (gtfFile, gtfFile), shell=True)
	subprocess.call("tabix -p gff -s 1 -b 4 -e 5 %s.gz" % gtfFile, shell=True)	# chromosome, start, end
maximumDistance = int(maximumDistance)

reader = openMAF(inputFile)
for paragraph in reader:
	alignments = []
	for line in paragraph.rstrip('\n').split('\n'):
		tokens = line.split()
		if tokens[0] == 's':
			(src, start, size, strand, srcSize, text) = tokens[1:]
			(start, size, srcSize) = map(int, (start, size, srcSize))
			if strand == '+':
				end = start + size
			if strand == '-':
				end = srcSize - start
				start = end - size
			# start, end: 0-based exclusive
			# Alignment: 1-based inclusive
			alignments.append(Alignment(src, start + 1, end, strand, text))
	for alignment in alignments:
		alignment.annotateGene(gtfFile)
	if all([alignment.index2genes for alignment in alignments]):
		length = min([len(alignment.sequence) for alignment in alignments])
		(startIndex, endIndex) = (0, 0)
		positionsList = []
		genesListList = []
		for index in range(length):
			bases = [alignment.sequence[index].upper() for alignment in alignments]
			if not '-' in bases and len(set(bases)) > 1:
				positions = [alignment.index2position[index] for alignment in alignments]
				genesList = [alignment.index2genes[index] if alignment.index2genes.has_key(index) else None for alignment in alignments]
				if all(genesList):
					if index >= endIndex:
						if startIndex < endIndex:
							printMarkers(alignments, startIndex, endIndex, positionsList, genesListList)
						startIndex = getStartIndex(alignments, index, max(index - maximumDistance, 0))
						positionsList = []
						genesListList = []
					endIndex = getEndIndex(alignments, index + 1, min(index + 1 + maximumDistance, length))
					positionsList.append(positions)
					genesListList.append(genesList)
		if startIndex < endIndex:
			printMarkers(alignments, startIndex, endIndex, positionsList, genesListList)
reader.close()
