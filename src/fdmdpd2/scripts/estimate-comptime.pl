#!/usr/bin/perl -w
# Average computation time per integration time step and per particle
# Arguments are directories!

use strict;

while(my $dir = shift) {
	my $min = 1e100;
	my $max =-1e100;
	my $maxstep = -1e100;
	my $minstep = 1e100;
	my $nN = 0;
	my $cpus = 0;

	next if (! -d $dir);
	
	opendir(DIR, $dir);
	while(my $file = readdir(DIR)) {
		my ($readtime, $writetime) = (stat("$dir/$file"))[8,9];
#print "$dir/$file: $readtime, $writetime\n";
		open(FH, "<$dir/$file") or die "Cannot open $dir/$file for reading: $!\n";
		my $line = <FH>;
		next if (!defined $line);
		if($line =~ /^# L=([\w\.]+) ([\w\.]+) ([\w\.]+) n=(\d+) N=(\d+) t=([\w\.]+)/) {
			$min = $writetime if ($writetime < $min || !defined($min));
			$max = $writetime if ($writetime > $max || !defined($max));
			$nN = $4 * $5;
#print "wt: $writetime $min $max\n";
		}
		close(FH);

		if($file =~ /output(\d+).dat/) {
			my $tmp = $1;
			$minstep = $tmp if($tmp < $minstep);
			$maxstep = $tmp if($tmp > $maxstep);
		}
	}
	closedir(DIR);

	my @list = <$dir.e*>;
	if (@list == 1) {
		open(FH, "<$list[0]") or die "Cannot open $list[0] for reading: $!\n";
		while(<FH>) {
			if(/^mpiexec: All (\d+) tasks \(spawn (\d+)\) started\./) {
				$cpus=$1;
			}
		}
		close(FH); 
	}
	
	my $diff=$max - $min;
	my $steps=$maxstep - $minstep;
	if($diff > 0 && $steps > 0 && $cpus > 0 && $nN > 0) {
		my $dt = $cpus * $diff / ($steps * $nN);
		print "$dir $dt $diff $nN $steps $cpus\n";
	}
}


# mpiexec: All 32 tasks (spawn 0) started.

# L=50 41.8827 41.8827 n=800 N=32 t=980.01
#
