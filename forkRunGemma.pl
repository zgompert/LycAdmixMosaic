#!/usr/bin/perl
#
# fit gemma BSLMM for Lycaeides wings 
#

use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);

$Nph = 17; ## 17 traits

foreach $ph (1..$Nph){ 
	foreach $ch (0..19){
		sleep 2;
		$pm->start and next;
		$g = "../geno/GNP.geno";
		$p = "../../traits/pheno_GNP.txt";
		$o = "o_lyc_GNP_ph$ph"."_ch$ch";
   		system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 500000 -s 1000000\n";
		$g = "../geno/SIN.geno";
		$p = "../../traits/pheno_SIN.txt";
		$o = "o_lyc_SIN_ph$ph"."_ch$ch";
   		system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 500000 -s 1000000\n";
		$g = "../geno/YBG.geno";
		$p = "../../traits/pheno_YBG.txt";
		$o = "o_lyc_YBG_ph$ph"."_ch$ch";
   		system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 500000 -s 1000000\n";
		$pm->finish;
	}
}
$pm->wait_all_children;

