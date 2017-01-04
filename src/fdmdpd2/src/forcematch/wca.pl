#!/usr/bin/perl -w
#
# WCA decomposition of an interaction potential
#
# Ref.: J. D. Weeks et al., J Chem Phys 54, 5237 (1971)
#
# The expected input format is "r U(r)\n"

use strict;

my $r_min = -1;
my $U_min = 1e10;
my @all_r;
my @all_U;

while(<>) {
	next if (/^#/);
	chomp;
	my ($r, $U) = split / /;
	if ($U < $U_min) {
		$r_min = $r;
		$U_min = $U;
	}
	push @all_r, $r;
	push @all_U, $U;
}

printf STDERR "r_min=$r_min, U_min=$U_min\n";

for (my $i = 0; $i < @all_r; $i++) {
	my $U_rep = ($all_r[$i] < $r_min) ? $all_U[$i] - $U_min : 0.;
	my $U_att = ($all_r[$i] < $r_min) ? $U_min : $all_U[$i];
	print "$all_r[$i] $U_rep $U_att\n";
}
