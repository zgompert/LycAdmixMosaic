#!/usr/bin/perl

## This script returns a genotype matirx (locus by ind) with genotype
## means (point estimates) from a genotype likelihood file; a HW prior is used based on allele frequencies provided in a separate file

## USAGE: perl gl2MaxGest.pl af_file.txt file.gl
use warnings;

$af = shift (@ARGV);
$in = shift (@ARGV);


## read in and store maf's
open (IN, $af) or die "read failed: $af\n";
while (<IN>){
	chomp;
	@line = split(/\s+/,$_);
	push (@af,$line[2]);
}
close (IN);


## read through gl file and estimate genotypes
open (IN, $in) or die "read failed: $in\n";
$out = $in;
$out =~ s/gl$/txt/;
open (OUT, "> pntest_$out") or die;
while (<IN>){
    chomp;
    if (s/^[0-9]+:\d+\s+//){ ## this line has genotype data, get rid of locus id
	$p = shift(@af); ## get alt. af. for this locus
        $prior[0] = 1 * ((1-$p) ** 2); ## for tets
        $prior[1] = 2 * $p * (1-$p) ;
        $prior[2] = 1 * ($p ** 2); 

	@line = split(" ",$_);
	@gest = ();
	while (@line){
	    $sum = 0;
	    for $i (0..2){ ## three genotyple likelihoods for each individual
		$gl[$i] = shift(@line);
		$gl[$i] = (10 ** ($gl[$i]/-10)) * $prior[$i];
		$sum += $gl[$i];
	    }
    	    for $i (0..2){ ## normalize
		$gl[$i] = $gl[$i]/$sum;
	    }
   	    $gest = 0;
	    $mgl = $gl[0]; ## set to 0
	    for $i (1..2){
		    if($gl[$i] > $mgl){
			    $mgl = $gl[$i];
			    $gest = $i;
		    }
	    }
		

    	    print OUT "$gest ";
	}
	print OUT "\n";	
    }
    else {
	print "failed to match $_\n";
    }
}
close (IN);
