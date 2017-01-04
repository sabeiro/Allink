#!/usr/bin/perl -w
#
# This program calculates the number composition fluctuations of particles in
# lateral boxes with a size of (L_y/n_y, L_z/n_z). n_y, n_z run from 1 to the
# specified value.
#
use strict;
use POSIX qw/floor/;

my @box;
my $nmax = shift or syntax();
my @L;
my $first = 1;
my @N;
my @N2;
my @cnt;
my $blocks;

while(<>) {
	chomp;
	if (/^# L=([^ ]+) ([^ ]+) ([^ ]+) t=([^ ]+) blocks=(\d+)/) {
		print STDERR "$_\n";

		compressibility() if ($first == 0);
		$first = 0;

		@box = ();
		
		$L[0] = $1;
		$L[1] = $2;
		$L[2] = $3;
		$blocks = $5;
		next;
	} elsif (/^#/) {
		next;
	}
	chomp;
	my @tmp = split / /;
	my $y = $tmp[1] - floor($tmp[1] / $L[1]) * $L[1];
	my $z = $tmp[2] - floor($tmp[2] / $L[2]) * $L[2];
	my $t = $tmp[6] % $blocks;

	$y = abs($y);
	$z = abs($z);

	die "y=$y, tmp[1]=$tmp[1] is invalid!\n\n" if ($y < 0. || $y > $L[1]);
	die "z=$z, tmp[2]=$tmp[2] is invalid!\n\n" if ($z < 0. || $z > $L[2]);
	die "t=$t is invalid!\n\n" if ($t < 0 || $t >= $blocks);
	
	my $a = floor($y / ($L[1] / $nmax));
	my $b = floor($z / ($L[2] / $nmax));
	$box[$a * $nmax + $b][$t]++;
}
compressibility();

#for (my $n = 1; $n < $nmax; $n++) {
#	my $by = $L[1] / $n;
#	my $bz = $L[2] / $n;
#	print "$by $bz $n";
#	for (my $k = 0; $k < $blocks / 2; $k++) {
#		if ($cnt[$n][$k] == 0) {
#			print "cnt=$cnt[$n][$k] for n=$n, k=$k, nmax=$nmax\n";
#		}
#		my $cnt = $cnt[$n][$k];
#		my $n1 = $N[$n][$k] / $cnt;
#		my $n2 = $N2[$n][$k] / $cnt;
#	
#		my $phi = $by * $bz * ($n2 / $n1 / $n1 - 1.);
#		print " $phi";
#	}
#	print "\n";
#}

sub compressibility
{
	for (my $i = 0; $i < $nmax; $i++) {
		for (my $j = 0; $j < $nmax; $j++) {
			for (my $k = 0; $k < $blocks; $k++) {
				my $tmp = $box[$i * $nmax + $j][$k];
				$tmp = 0 if (!defined $tmp);
				print "$tmp ";
			}
			print "\n";
		}
	}
}

sub syntax
{
	print STDERR "Syntax: $0 [nmax]\n\n";
	exit 0;
}


