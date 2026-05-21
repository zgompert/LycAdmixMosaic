#!/usr/bin/perl
#
# correct the chromosome ids, iSMC numbers them sequentially,

# get order of chromosomes

open(IN,"chOrd.txt") or die;
while(<IN>){
	chomp;
	s/Ch//;
	push(@ch,$_);
}
close(IN);

## give bedgraph files
## idas_ismc.rho.250kb.bedgraph

foreach $bed (@ARGV){
	$out = "ch_$bed";
	open(IN,$bed) or die;
	open(OUT, "> $out") or die;
	while(<IN>){
		chomp;
		if(m/^chr(\d+)/){
			$ach = $ch[$1-1];
			s/^chr\d+/$ach/ or die "failed here: $_ :: $ach\n";
		}
		print OUT "$_\n";
	}
	close(IN);
	close(OUT);
}

