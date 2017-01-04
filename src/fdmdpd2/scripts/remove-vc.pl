#!/usr/bin/perl -w
# Removes center of mass velocity for each block separately
#
#
use strict;

my $nN;
my @coord = ();
my $i = 0;
my $block = 0;

while (<>) {
	if (/^# n=(\d+) N=(\d+) name=(\w+)/) {
		$nN = $1 * $2;
		$block++;
		print;
	} elsif (/^# /) {
		print;
	} else {
		chomp;
		my @tmp = split / /;
		for (my $j = 0; $j < @tmp; $j++) {
			$coord[$i][$j] = $tmp[$j];
		}
		$i++;
	}

	if ($i > 0 && $i == $nN) {
		my @v = ();

		for ($i = 0; $i < $nN; $i++) {
			$v[0] += $coord[$i][3];
			$v[1] += $coord[$i][4];
			$v[2] += $coord[$i][5];
		}

		$v[0] /= $nN;
		$v[1] /= $nN;
		$v[2] /= $nN;

		print STDERR "v_$block = (", join(' ', @v), ")\n";

		for ($i = 0; $i < $nN; $i++) {
			$coord[$i][3] -= $v[0];
			$coord[$i][4] -= $v[1];
			$coord[$i][5] -= $v[2];
		}

		for ($i = 0; $i < $nN; $i++) {
			print join(' ', @{$coord[$i]}), "\n";
		}

		@coord = ();
		$i = 0;
	}
}

