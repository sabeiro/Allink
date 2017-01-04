#!/usr/bin/perl -w
use strict;
my $m = shift or die "Mean value missing";
my $a = 0.;
my $a2 = 0.;
my $n = 0;

while(<>) {
	my @tmp = split / /;
	$a2 += $tmp[1] * $tmp[1];
	$n++;
}

$a2 = sqrt($a2 / $n - $m * $m);
$a2 /= sqrt($n);
print "$a2\n";
