# SEXCMD : SEX Created Marker and Determination

 USAGE : perl determine.pl sex_marker_filtered.hg38.final.fasta input.fastq.gz 2> /dev/null 1> input.SEX_OUTPUT

 SEXCMD is available at https://github.com/lovemun/SEXCMD

 The code is written in Perl and bash script. This tool is supported on Linux and needs Python, lastz, blastall and bwa. The tool uses gzip compressed fastq file as input and R1 and R2 file can be used in case of paired-end or mate pair. And you should give marker fasta file path, input fastq file path as auguements.

# Required
 1. Linux Program
 2. Python
 3. bwa(version 7.0 or later)
 4. blastall
 5. lastz(Release 1.02.00)

# Example
 If you want to create marker file, you can make own sex marker file for your datafiles. There are 6 process steps for creating own sex marker file.

1. Extract sex chromosome from Reference genome(ex. Human)

  python util/0.extract_chr_from_REF.py hg38.fasta X > hg38.chrX.fasta

  python util/0.extract_chr_from_REF.py hg38.fasta Y > hg38.chrY.fasta


2. Mapping between sex chromosomes

 lastz hg38.chrX.fasta hg38.chrY.fasta --format=maf > chrX_chrY.maf

3. Generate sex marker text file

 python util/1.sex_marker.py chrX_chrY.maf hg38.ucsc.gtf 35 > sex_marker.hg38.txt 
 (distance between marker : 35(recommended))

 python util/2.sex_marker_filtered.py sex_marker.hg38.txt hg38.ucsc.gtf 151 5 > sex_marker_filtered.hg38.txt 
 (read length : 151 admit missmatch count : 5)

4. Convert sex marker file format(text -> fasta)

 python util/3.make_markerFASTA.py sex_marker_filtered.hg38.txt

5. Filtered ONLY MAPPED to chrX, chrY

 blastn -db hg38.fa -query sex_marker_filtered.hg38.fasta -num_threads 4 -word_size 8 -outfmt 7 > sex_marker_filtered.hg38.fasta.blastm9
 
 python util/4.blast_filter.py sex_marker_filtered.hg38.fasta sex_marker_filtered.hg38.fasta.blastm9 sex_marker_filtered.hg38.final.fasta

6. Index sex marker file

 bwa index -a is sex_marker_filtered.hg38.final.fasta

Final sex marker : sex_marker_filtered.hg38.final.fasta
