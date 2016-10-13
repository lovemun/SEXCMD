# Author: Jiwoong Kim (jiwoongbio@gmail.com) supervised by Namshin Kim (n@rna.kr)
# type : XY or ZW
use strict;
use warnings;

my ($fastaFile,$type,@fastqFileList) = @ARGV;
if(@fastqFileList){
	chomp(my @lineList = map {`bwa aln -t 4 $fastaFile $_ | bwa samse $fastaFile - $_ | samtools view -S -F 4 -q 30 - | cut -f3 | sort | uniq -c | awk '{print \$2"\\t"\$1}'`} @fastqFileList);
	my %readCountHash = map {$_->[0] => $_->[1]} map {[split(/\t/, $_)]} @lineList;
	my %genderCountHash = ();
	open(my $reader, "grep '^>' $fastaFile | sed 's/^>//' | paste - - |");
	while(my $line = <$reader>) {
		chomp($line);
		my ($readCountX, $readCountY) = map {defined($_) ? $_ : 0} @readCountHash{my ($nameX, $nameY) = split(/\t/, $line)};
		if($type eq "XY") {
			$genderCountHash{my $gender = getGender_XY($readCountX, $readCountY)}++;
			print join("\t", $nameX, $nameY, $readCountX, $readCountY, $gender), "\n";
		} elsif($type eq "ZW") {
			$genderCountHash{my $gender = getGender_ZW($readCountX, $readCountY)}++;
			print join("\t", $nameX, $nameY, $readCountX, $readCountY, $gender), "\n";
		} else {
			print STDERR			""
		}
	}
	close($reader);
	my ($genderCountF, $genderCountM) = map {defined($_) ? $_ : 0} @genderCountHash{'F', 'M'};
	if($genderCountF < $genderCountM) {
		print "M\n";
	} elsif($genderCountF > $genderCountM) {
		print "F\n";
	} else {
		print "?\n";
	}
} else {
	print STDERR '
Usage: perl determineGender.pl [genderMarker.fasta] [sex_system(XY, ZW, ..)] [input.fastq ...]
       - The command "bwa" should be available in the current working directory.
       - genderMarker.fasta should be indexed by "bwa".

';
}

sub getGender_XY {
	my ($readCountX, $readCountY) = @_;
	if($readCountX + $readCountY > 0) {
		if($readCountY > 0) {
			return 'M';
		} else {
			return 'F';
		}
	} else {
		return '?';
	}
}

sub getGender_ZW {
	my ($readCountZ, $readCountW) = @_;
	if($readCountZ + $readCountW > 0) {
		if($readCountW > 0) {
			return 'F';
		} else {
			return 'M';
		}
	} else {
		return '?';
	}
}
