#!/usr/bin/perl
#
# only keep chromosomes, and rename them, header and SNPs
#

## get chroms
open(IN,"lgs.txt");
while(<IN>){
	chomp;
	@line = split(/\s+/,$_);
	$chrom{$line[0]} = $line[1];
}
close(IN);

foreach $vcf (@ARGV){
	$out = "clean_$vcf";
	open(IN, $vcf);
	open(OUT, "> $out");
	while(<IN>){
		chomp;
		if(m/\#/){
			if(m/Scaffold_(\d+)/){
				$sc = $1;
				if(defined($chrom{$sc})){
					s/Scaffold_\d+;HRSCAF_\d+/Ch$chrom{$sc}/;
					print OUT "$_\n";
				}
			} else{
				print OUT "$_\n";
			} 
		} elsif(m/^Scaffold_(\d+)/){
			$sc = $1;
                	if(defined($chrom{$sc})){
                        	s/Scaffold_\d+;HRSCAF_\d+/Ch$chrom{$sc}/;
                                print OUT "$_\n";
                        }
		}
	}
	close(IN);
	close(OUT);
}
