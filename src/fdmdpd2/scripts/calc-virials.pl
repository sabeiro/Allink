#!/usr/bin/perl -w
use strict;

my $rho = shift or syntax();
my $kN = shift or syntax();
my $chiN = shift or syntax();

my $v_aa = -2. * ($kN + 3.) / $rho;
my $v_bb = 0.1;
my $v_ab = $chiN / $rho + .5 * ($v_aa + $v_bb);
my $w_aaa = 1.5 * ($kN + 2.) / ($rho * $rho);
my $w_aab = $w_aaa;
my $w_abb = $w_aaa;
my $w_bbb = 0.0;
printf("# v=%lg %lg %lg w=%lg %lg %lg %lg\n",
					$v_aa, $v_ab, $v_bb,
					$w_aaa, $w_aab, $w_abb, $w_bbb);

sub syntax
{
	print "Syntax: calc-virials.pl [rho_coex] [kappaN] [chiN]\n\n";
	exit(0);
}
