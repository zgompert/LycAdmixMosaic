#!/usr/bin/perl
#
## version of treemix fork script that obtains the ML result across $N runs


use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);

$N = 20;

foreach $fi (@ARGV){
	$fi =~ m/(\d+)/ or die "failed here $fi\n";
	$ch = $1;
	foreach $m (0..8){
		$pm->start and next;
		$out = "tro_ch$ch"."_m$m";
		system "treemix -i $fi -o $out -k 100 -m $m -root MEN12\n";
		open(IN, "$out\.llik");
		$a = <IN>;
		$a = <IN>;
		chomp($a);
		$a =~ m/:\s+(\-[0-9\.]+)/ or die "can't find the ll: $a\n";
		$ll = $1;
		close(IN);
		print "starting ll = $ll\n";
		foreach $i (1..$N){
			$tout = "temp_tro_ch$ch"."_m$m";
			system "treemix -i $fi -o $tout -k 100 -m $m -root MEN12\n";
			open(IN, "$tout\.llik");
			$a = <IN>;
			$a = <IN>;
			chomp($a);
			$a =~ m/:\s+(\-[0-9\.]+)/ or die "can't find the ll: $a\n";
			$llt = $1;
			if($llt > $ll){
				$ll = $llt;
				system "mv $tout\.treeout.gz $out\.treeout.gz\n";
				system "mv $tout\.modelcov.gz $out\.modelcov.gz\n";
				system "mv $tout\.vertices.gz $out\.vertices.gz\n";
				system "mv $tout\.cov.gz $out\.cov.gz\n";
				system "mv $tout\.covse.gz $out\.covse.gz\n";
				system "mv $tout\.edges.gz $out\.edges.gz\n";
				system "mv $tout\.llik $out\.llik\n";
			}
		}
		$pm->finish;
	}
}

$pm->wait_all_children;


