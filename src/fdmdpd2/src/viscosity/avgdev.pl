#!/usr/bin/perl -w
use strict;
my $m = shift or die "Mean value missing";
my $a = 0.;
my $n = 0;

while(<>) {
	my @tmp = split / /;
	$a += abs($tmp[1] - $m);
	$n++;
}
$a /= $n;
print "$a\n";
