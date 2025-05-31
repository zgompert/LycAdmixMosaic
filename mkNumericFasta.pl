#!/usr/bin/perl

## this is to figure out which SNPs are variable in the fasta

foreach $i (1..23){
	system "grep -v \"^>\" sub_max_chrom$i.fasta | perl -p -i -e 'tr/ACGTN/12345/' | sed 's/./& /g' > text_max_chrom$i.fasta\n";
	#system "grep -v \"^>\" sub_chrom$i.fasta | perl -p -i -e 'tr/ACGT/1234/' | sed 's/./& /g' > text_chrom$i.fasta\n";

}
