#!/usr/bin/perl

## provide snp file and  param files
$snpfile = shift(@ARGV);
open(SNP, $snpfile) or die "could not open the SNP file\n";
while(<SNP>){
	chomp;
	push(@snps,$_);
}
close(SNP);


foreach $in (@ARGV){
	foreach $snp (@snps){ ## define and 0 all SNP effects
		$eff{$snp} = 0;
	}
## combine information over reps and sort
	open(OUT, "> pip_$in") or die "failed to write mav $in\n";
	$nreps=20;
	foreach $rep (0..19){
		$in =~ s/ch\d+/ch$rep/ or die "failed sub for $in to $rep\n";
		print "working on $in\n";
		open(IN, $in) or die "failed to open the infile $in\n";
		<IN>; ## burn header
		while(<IN>){
			chomp;
			@line = split(/\s+/,$_);
			$eff{$line[1]} +=  $line[6] /$nreps;
		}
		close(IN);
	}

	foreach $snp (@snps){
		print OUT "$snp $eff{$snp}\n";
	}
	close(OUT);
}
