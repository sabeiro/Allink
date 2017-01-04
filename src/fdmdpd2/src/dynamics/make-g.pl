#!/usr/bin/perl -w
use strict;

my @dp = ();
my $c = 0;
my $dq = 0.05;
my @histo = ();
my @histo2 = ();
my @histoc = ();

while(<>) {
	next if (/^\D/);
	# q1 q2 |q| S
	my ($qy, $qz, $q, $S) = split / /;
	my $bin = int($q / $dq);
	$histo[$bin] += $S;
	$histo2[$bin] += $S * $S;
	$histoc[$bin]++;
}

for (my $i = 0; $i < @histo; $i++) {
	my $q = $i * $dq;
	next if (!(defined $histo[$i]));
	my $h = $histo[$i] / $histoc[$i];
	my $h2 = $histo2[$i] / $histoc[$i];
	my $n = $histoc[$i];

	print "$q $h $h2 $n\n";
}
