#!/usr/bin/perl -w
#
# This scripts performs a "Galilei-boost" into the reference frame of the
# center-of-mass of the given system.
#
# Syntax: gboost.pl [input.dat] > [output.dat]
#
use strict;

my $fn = shift;
my @vc;
my $count = 0;

open(FH, "<$fn") or die "Cannot open file $fn for reading: $!\n";

while(<FH>) {
	next if (/^#/);
	my @tmp = split(/ /);
	$vc[0] += $tmp[3];
	$vc[1] += $tmp[4];
	$vc[2] += $tmp[5];
	$count++;
}

$vc[0] /= $count;
$vc[1] /= $count;
$vc[2] /= $count;
print STDERR "vc=". join(' ', @vc) ."\n";

seek(FH, 0, 0);
while(<FH>) {
	if (/^#/) {
		print;
	} else {
		my @tmp = split(/ /);
		$tmp[3] -= $vc[0];
		$tmp[4] -= $vc[1];
		$tmp[5] -= $vc[2];
		print join(' ', @tmp);
	}
}
close(FH);
