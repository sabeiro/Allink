#!/usr/bin/perl -w
use strict;

my $factor = shift or die "Stretching factor missing\n";

while (<>) {
	if (/^# L=([^ ]+) ([^ ]+) ([^ ]+) t=([^ ]+) blocks=(\d+)/) {
		my $Ly = $2 * $factor;
		my $Lz = $3 * $factor;
		print "# L=$1 $Ly $Lz t=$4 blocks=$5\n";
	} elsif (/^#/) {
		print;
	} else {
		my @tmp = split / /;
		$tmp[1] *= $factor;
		$tmp[2] *= $factor;
		print join(' ', @tmp);
	}
}
