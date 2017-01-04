#!/usr/bin/perl -w
# 
# Split a single output file containing many snapshots into many small output
# files, that contain only one snapshot.
#
# Syntax: ./splitoutput.pl [output.dat]
# 
use strict;

my $count = 0;

# L=50 71.9942 71.9942 t=0.00000000 blocks=1
# v=-11.8571 -3.73571 0.1 w=0.627551 0.627551 0.627551 0
# a2=0.9 a3=1 Re=3.5 N=16 ks=19 kb=5 l0=0
# n=18720 N=16 name=LIPID

while(my $line = <>) {
	my $fn = sprintf("split%05d.dat", $count++);
	open(FH, ">$fn");

	$line =~ /blocks=(\d+)/;
	my $blocks = $1;
	print FH $line;
	
	$line = <>;
	print FH $line;
	
	$line = <>;
	print FH $line;

	for (my $i = 0; $i < $blocks; $i++) {
		$line = <>;
		$line =~ /n=(\d+) N=(\d+)/;
		my $n = $1;
		my $N = $2;
		print FH $line;
		for (my $j = 0; $j < $n * $N; $j++) {
			$line = <>;
			print FH $line;
		}
	}

	close(FH);
}
