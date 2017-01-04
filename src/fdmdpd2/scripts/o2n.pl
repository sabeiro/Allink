#!/usr/bin/perl -w
use strict;

my $n;
my $N;

while(<>) {
	if (/^# L=/) {
		/# L=([\d\.]+) ([\d\.]+) ([\d\.]+) n=(\d+) N=(\d+) t=([\d\.]+)/;

		my $lx = $1;
		my $ly = $2;
		my $lz = $3;
		$n = $4;
		$N = $5;
		my $t = $6;

		print "# L=$lx $ly $lz t=$t blocks=1\n";
	} elsif (/^# v=/) {
		print;
		print "# a2=0.9 a3=1.0 Re=3.5 N=16 ks=19.0 kb=5.0 l0=0\n";
		print "# n=$n N=$N name=LIPID\n";
	} else {
		print;
	}	
}
