#!/usr/bin/perl
#
# run lepmap3 order module
#

foreach $ch (1..23){
	print "Working on $ch\n";
	$ord = "chr_$ch.physical_order.txt";
	$map = "chr_$ch.map";
	$log = "chr_$ch.log";
	#system "java -cp ~/source/lm3/bin/ OrderMarkers2 data=f_lmdata.call evaluateOrder=$ord improveOrder=0 recombination2=0 outputPhasedData=1 minError=0.01 proximityScale=100 \"scale=0.5M/N\" > $map 2> $log\n";
	#system "java -cp ~/source/lm3/bin/ OrderMarkers2 data=f_lmdata.call evaluateOrder=$ord improveOrder=0 recombination2=0 outputPhasedData=1 minError=0.01 proximityScale=100 \"scale=0.5M/N 2\" > $map 2> $log\n";
	system "java -cp ~/source/lm3/bin/ OrderMarkers2 data=f_lmdata.call evaluateOrder=$ord improveOrder=0 recombination2=0 outputPhasedData=1 minError=0.1 proximityScale=100 \"scale=0.5M/N 2\" > $map 2> $log\n";

}

