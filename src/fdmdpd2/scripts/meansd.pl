#!/usr/bin/perl -w
use strict;

my $n = 0;
my $m = ();
my $m2 = ();

while(<>) {
	chop;
	next if (/^#/);
	my @tmp = split / /;
	$m += $tmp[0];
	$m2 += $tmp[0] * $tmp[0];
	$n++;
}
$m /= $n;
my $sd = sqrt($m2/$n - $m * $m);
my $sdm = $sd / sqrt($n - 1.);
print "$m $sd $sdm $n\n";
