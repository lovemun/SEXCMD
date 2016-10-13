#!/usr/bin/env python

import sys
import os.path
import subprocess

class Marker:
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
			indexList = range(len(self.sequence))[::-1]
		position = self.start
		for index in indexList:
			if not self.sequence[index] == '-':
				self.index2position[index] = position
				self.position2index[position] = index
				position += 1
		self.mismatchPositions = []
		self.genes = set()
	# Annotate genes using tabix
	def annotateGene(self, gtfFile):
		for mismatchPosition in self.mismatchPositions:
			pipe = subprocess.Popen("tabix %s.gz %s:%d-%d" % (gtfFile, self.chromosome, mismatchPosition, mismatchPosition), shell=True, stdout=subprocess.PIPE).stdout
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
							self.genes.add(value)
			pipe.close()
	def printLine(self):
		print '\t'.join([self.chromosome, str(self.start), str(self.end), self.strand, self.sequence, ','.join(map(str, self.mismatchPositions)), ','.join(self.genes)])

def getMarker(line):
	tokens = line.rstrip('\n').split('\t')
	marker = Marker(chromosome = tokens[0], start = int(tokens[1]), end = int(tokens[2]), strand = tokens[3], sequence = tokens[4])
	marker.mismatchPositions = map(int, tokens[5].split(','))
	marker.genes = set(tokens[6].split(','))
	return marker

def getMarkerPart(marker, startIndex, endIndex):
	(start, end) = sorted([marker.index2position[startIndex], marker.index2position[endIndex - 1]])
	markerPart = Marker(marker.chromosome, start, end, marker.strand, marker.sequence[startIndex:endIndex])
	markerPart.mismatchPositions = [mismatchPosition for mismatchPosition in marker.mismatchPositions if start <= mismatchPosition <= end]
	markerPart.genes = marker.genes
	return markerPart

def stripMarkerList(markerList):
	length = min([len(marker.sequence) for marker in markerList])
	startIndex = 0
	while startIndex < length:
		if len(set([marker.sequence[startIndex] for marker in markerList])) == 1:
			break
		startIndex += 1
	endIndex = length
	while endIndex > 0:
		if len(set([marker.sequence[endIndex - 1] for marker in markerList])) == 1:
			break
		endIndex -= 1
	if startIndex < endIndex:
		return [getMarkerPart(marker, startIndex, endIndex) for marker in markerList]
	else:
		return []

# Input
(inputFile, gtfFile, minimumLength, minimumMismatch) = sys.argv[1:]
if not os.path.isfile("%s.gz.tbi" % gtfFile):
	subprocess.call("grep -v '^#' %s | sort --field-separator='\t' -k1,1 -k4,4n -k5,5n | bgzip > %s.gz" % (gtfFile, gtfFile), shell=True)
	subprocess.call("tabix -p gff -s 1 -b 4 -e 5 %s.gz" % gtfFile, shell=True)
(minimumLength, minimumMismatch) = map(int, [minimumLength, minimumMismatch])

reader = open(inputFile, 'r')
chromosomePositionCountDict = {}
for line in reader:
	marker = getMarker(line)
	for mismatchPosition in marker.mismatchPositions:
		chromosomePositionCountDict.setdefault(marker.chromosome, {}).setdefault(mismatchPosition, 0)
		chromosomePositionCountDict[marker.chromosome][mismatchPosition] += 1
reader.seek(0)
lines = filter(lambda x: x, [reader.readline(), reader.readline()])
while lines:
	markerList = [getMarker(line) for line in lines]
	(indexList, ) = set([tuple([marker.position2index[mismatchPosition] for mismatchPosition in marker.mismatchPositions]) for marker in markerList])
	startIndex = 0
	for index in indexList:
		if any([chromosomePositionCountDict[marker.chromosome][marker.index2position[index]] > 1 for marker in markerList]):
			endIndex = index
			if endIndex - startIndex >= minimumLength:
				for marker in stripMarkerList([getMarkerPart(marker, startIndex, endIndex) for marker in markerList]):
					if marker is not None and len(marker.sequence) >= minimumLength and len(marker.mismatchPositions) >= minimumMismatch:
						marker.annotateGene(gtfFile)
						marker.printLine()
			startIndex = index + 1
	endIndex = len(marker.sequence)
	if endIndex - startIndex >= minimumLength:
		for marker in stripMarkerList([getMarkerPart(marker, startIndex, endIndex) for marker in markerList]):
			if marker is not None and len(marker.sequence) >= minimumLength and len(marker.mismatchPositions) >= minimumMismatch:
				marker.annotateGene(gtfFile)
				marker.printLine()
	lines = filter(lambda x: x, [reader.readline(), reader.readline()])
reader.close()
