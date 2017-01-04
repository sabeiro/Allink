#!/usr/bin/perl -w
use strict;

my $N;
my @cm = ();
my $na = 0;
my $nb = 0;

while(<>) {
	if (/^# n=(\d+) N=(\d+) name=(\w+)/) {
		$N = $2;
		print "# n=$1 N=1 name=$3\n";
	} elsif (/^#/) {
		print;
	} else {
		chop;
		my @tmp = split / /;
		if ($tmp[6] == 0) {
			for my $i (0 .. 5) {
				$cm[$i] += $tmp[$i];
			}
			$na++;
		} else {
			$nb++;
		}

		if ($na + $nb == $N) {
			for my $i (0 .. 5) {
				$cm[$i] /= $na;
			}
			printf('%lg ', $_) foreach (@cm);
			print "0\n";
			@cm = ();
			$na = 0;
			$nb = 0;
		}
	}
}

