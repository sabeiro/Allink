#!/usr/bin/perl -w
use strict;

my @r;
my @ur;
my @uz;
my $c = 0;

while(<>) {
	next if(/^#/);
	my @tmp = split /\s/;
	$r[$c] = $tmp[0];
	$ur[$c] = $tmp[1];
	$uz[$c] = $tmp[3];
	$c++;
}

my $zmin = -20.; #-1. * $r[$c-1];
my $zmax = 20.; #$r[$c-1];
#my $dr = $r[1] - $r[0];
my $dr = 0.01;

for (my $i = 0; $i < $c; $i++) {
	for (my $j = $c - 1; $j > 0; $j--) {
		my $z = -$r[$j];
		my $r = $r[$i];
		my $u = $ur[$i] * $r / sqrt($r * $r + $z * $z);
		$u += -1. * $uz[$j] * $z / sqrt($r * $r + $z * $z);
		printf("%lg\t%lg\t%lg\n", $r, $z, $u);
	}
	for (my $j = 0; $j < $c; $j++) {
		my $z = $r[$j];
		my $r = $r[$i];
		my $u;
		if ($i == 0 && $j == 0) {
			$u = 0.;
		} else {
			$u = $ur[$i] * $r / sqrt($r * $r + $z * $z);
			$u += $uz[$j] * $z / sqrt($r * $r + $z * $z);
		}
		printf("%lg\t%lg\t%lg\n", $r, $z, $u);
	}
	print "\n";
}

