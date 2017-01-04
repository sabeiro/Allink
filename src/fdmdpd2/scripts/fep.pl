#!/usr/bin/perl -w
use strict;

my @du = ();
my $count = 0;

while(<>) {
	next if(!/^FEP: /);
	s/^FEP: //;
	s/\n$//;
	my @tmp = split / /;
	for (my $i = 0; $i < @tmp; $i++) {
		my $x = $tmp[$i];
		$du[$i] .= "$x ";
	}
	$count++;
}

open(PARI, "|gp") or die "Cannot open Pari calculator: $!\n";

for (my $i = 0; $i < @du; $i++) {

	my @tmp = split / /, $du[$i];
	my $mean = 0.;

	foreach my $x (@tmp) {
		$mean += $x;
	}
	$mean /= @tmp;
	#print "mean=$mean i=$i\n";
	
	my $s = '-1. * log((0';
	foreach my $x (@tmp) {
		$s .= "+exp(-1.*($x - $mean))";
	}
	$s .= ")/$count) + $mean\n";
	print PARI $s;
}

close(PARI);

