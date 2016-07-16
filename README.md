SEXCMD : SEX Created Marker and Determination

USAGE : perl determine.pl sex_marker_filtered.hg19.final.fasta input.fastq.gz

SEXCMD is available at https://github.com/pkoiu87/SEXCMD

The code is written in Perl and bash script. This tool is supported on Linux and needs Python, lastz, blastall and bwa. The tool uses gzip compressed fastq file as input and R1 and R2 file can be used in case of paired-end or mate pair. And you should give marker fasta file path, input fastq file path as auguements.

Required
1.Linux Program
2.Python
3.bwa(version 7.0 or later)
4.blastall
5.lastz(Release 1.02.00)

If you want to create marker file, you can make own sex marker file for your datafiles. There are 6 process steps for creating own sex marker file.
1.Extract sex chromosome from Reference genome(ex. Human)

python util/0.extract_chr_from_REF.py hg19.fasta X > hg19.chrX.fasta python util/0.extract_chr_from_REF.py hg19.fasta Y > hg19.chrY.fasta
1.Mapping between sex chromosomes

lastz hg19.chrX.fasta hg19.chrY.fasta --format=maf > chrX_chrY.maf
1.Generate sex marker text file

python util/1.sex_marker.py chrX_chrY.maf HG19.ucsc.gtf 35 > sex_marker.hg19.txt distance between marker : 35(recommended) python util/2.sex_marker_filtered.py sex_marker.hg19.txt HG19.ucsc.gtf 151 5 > sex_marker_filtered.hg19.txt read length : 151 admit missmatch count : 5
1.Convert sex marker file format(text -> fasta)

python util/3.make_markerFASTA.py sex_marker_filtered.hg19.txt
1.Filtered ONLY MAPPED to chrX, chrY

blastall -p blastn -d hg19 -i sex_marker_filtered.hg19.fasta -a 16 -W 8 -m 9 > sex_marker_filtered.hg19.fasta.blastm9 python util/4.blast_filter.py sex_marker_filtered.hg19.fasta sex_marker_filtered.hg19.fasta.blastm9 sex_marker_filtered.hg19.final.fasta
1.index sex marker file

bwa index -a is sex_marker_filtered.hg19.final.fasta

final sex marker : sex_marker_filtered.hg19.final.fasta
