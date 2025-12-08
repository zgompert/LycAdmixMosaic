#!/usr/bin/perl
#
# conver sam to bam, then sort and index 
#


use Parallel::ForkManager;
my $max = 16;
my $pm = Parallel::ForkManager->new($max);

foreach $iter (1..3){
	$pm->start and next; ## fork
        $out = "po_$iter.hdf5";
	system "popmod -i filtered_tcris_fha2013_bin_gs.gl -n 6000 -b 1000 -t 2 -o $out\n";
	$pm->finish;
}

$pm->wait_all_children;



