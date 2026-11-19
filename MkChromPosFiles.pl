#!/usr/bin/perl
#
# make chromosome specific input files
#

## first get scaffold to chromosome mappings 
open(ID, "LgScaffoldNames.txt") or die "Could not open the mapping file\n";
while(<ID>){
	chomp;
	m/(\S+)\s+(\S+)/ or die "failed here ID = $_\n";
	$chid{$1} = $2;
}
close(ID);

## open one file per chromosome for order and position
foreach $ch (sort keys %chid){
	$fh = "F"."$ch";
	$fh2 = "M"."$ch";
	$chn = $chid{$ch};
	$out = "chr_$chn".".physical_order.txt";
	print "$out\n";
	open ($fh, "> $out") or die "Could not write $id\n";
	$files{$ch} = $fh;

	$out = "chr_$chn".".pos";
	print "$out\n";
	open ($fh2, "> $out") or die "Could not write $id\n";
	$files2{$ch} = $fh2;
}

## now parse positions
$mcnt = 0; 
open(IN,"f_lmdata.call") or die;
while(<IN>){
	chomp;
	if(m/^(Sc\S+)\s+(\S+)/){
		$pos = $2;
		$chi = $1;
		$ch = $chid{$chi};
		print {$files2{$chi}} "$pos\n";
		$mcnt++;
		print {$files{$chi}} "$mcnt\n";
	}
}
close(IN);

## close the files
foreach $ch (sort keys %chid){
	close($files{$ch});
	close($files2{$ch});
}
