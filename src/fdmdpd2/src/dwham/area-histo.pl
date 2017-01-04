#!/usr/bin/perl -w
#
# Simple program to calculate the mean area as a function of the unnormalized
# order parameter without reweighting.
#
# Syntax: ./area-histo.pl [logfile 1] [logfile 2] ...
#

use strict;

my $bw = 0.002;
my @count;
my @histo_A;
my @histo_A2;

while (my $fn = shift)
{
	print STDERR "Opening '$fn'... \n";
	open(FH, "<$fn") or die "Cannot open '%s' for reading: $!\n";

	my $A;

	while(<FH>) {
		# get area
#		if (/sum=([^ ]+)/) {
#			$A = $1;
#		}
		if (/^l=([^ ]+) ([^ ]+) ([^ ]+)/) {
			$A = $2 * $3;
		}
		if (/^Sxx=([^ ]+)/) {
			my $bin = int($1 / $bw + .5);
			$histo_A[$bin] += $A;
			$histo_A2[$bin] += $A * $A;
			$count[$bin]++;
		}
	}
	close(FH);
}

print STDOUT "# S A A^2\n";
for (my $bin = 0; $bin * $bw < 1.; $bin++) {
	if (defined $count[$bin]) {
		my $S = ($bin + .5) * $bw;
		my $A = $histo_A[$bin] / $count[$bin];
		my $A2 = $histo_A2[$bin] / $count[$bin];
		print STDOUT "$S $A $A2\n";
	}
}


# time=9999.9 inter=-273958 intra=162849 kin=113116 sum=2005.99 <kT>=1.00708

