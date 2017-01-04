#!/usr/bin/perl -w
#
# This program calculates the number fluctuations of particles in lateral boxes
# with a size of (L_y/n_y, L_z/n_z). n_y, n_z run from 1 to the specified value.
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
		next;
	} elsif (/^#/) {
		next;
	}
	my @tmp = split / /;
	my $y = abs($tmp[1] - floor($tmp[1] / $L[1]) * $L[1]);
	my $z = abs($tmp[2] - floor($tmp[2] / $L[2]) * $L[2]);

#	die "y=$y, tmp[1]=$tmp[1] is invalid!\n\n" if ($y < 0. || $y > $L[1]);
#	die "z=$z, tmp[2]=$tmp[2] is invalid!\n\n" if ($z < 0. || $z > $L[2]);
	
	for (my $n = 1; $n < $nmax; $n++) {
		my $a = floor($y / ($L[1] / $n));
		my $b = floor($z / ($L[2] / $n));
		$box[$n][$a * $n + $b]++;
	}
}
compressibility();

for (my $n = 1; $n < $nmax; $n++) {
	if ($cnt[$n] == 0) {
		print "cnt=$cnt[$n] for n=$n, nmax=$nmax\n";
	}
	my $n1 = $N[$n] / $cnt[$n];
	my $n2 = $N2[$n] / $cnt[$n];
	my $cnt = $cnt[$n];

	my $by = $L[1] / $n;
	my $bz = $L[2] / $n;
	my $k = $by * $bz * ($n2 / $n1 / $n1 - 1.);
	print join(' ', $by, $bz, $n, $n1, $n2, $cnt, $k), "\n";
}

sub compressibility
{
	for (my $n = 1; $n < $nmax; $n++) {
		for (my $i = 0; $i < $n; $i++) {
			for (my $j = 0; $j < $n; $j++) {
				my $tmp = $box[$n][$i * $n + $j];
				$tmp = 0 if (!defined $tmp);
				$N[$n] += $tmp;
				$N2[$n] += $tmp * $tmp;
				$cnt[$n]++;
			}
		}
	}
}

sub syntax
{
	print STDERR "Syntax: $0 [nmax]\n\n";
	exit 0;
}


