#!/usr/bin/perl -w
use strict;

my @kN = qw/100./;
my @chiN = qw/30./;
my @header = ();
my @rho = ();

for (my $i = 15.0; $i <= 19.0; $i+=0.02) {
	push @rho, $i;
}

foreach my $kN (@kN) {
	foreach my $chiN (@chiN) {
		foreach my $rho (@rho) {
			my $v_aa = -2. * ($kN + 3.) / $rho;
			my $v_bb = 0.1;
			my $v_ab = $chiN / $rho + .5 * ($v_aa + $v_bb);
			my $w_aaa = 1.5 * ($kN + 2.) / ($rho * $rho);
			my $w_aab = $w_aaa;
			my $w_abb = $w_aaa;
			my $w_bbb = 0.0;
			my $header = sprintf("v=%lg %lg %lg w=%lg %lg %lg %lg\n",
					$v_aa, $v_ab, $v_bb,
					$w_aaa, $w_aab, $w_abb, $w_bbb);
			push @header, $header;
		}
	}
}

print "# ". @header ."\n";
print foreach (@header)
