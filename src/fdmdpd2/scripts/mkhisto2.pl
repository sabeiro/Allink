#!/usr/bin/perl -w
#
# Perl script to make a histogram with error bars
#
# Syntax: ./mkhisto.pl [bw] {file}
#
use strict;

my %histo = ();
my %histo2 = ();
my %histoc = ();

while(<>) {
	next if (/^#/);
	my @tmp = split / /;
	my $x = $tmp[0];
	my $y = $tmp[1];

	$histo{$x} += $y;
	$histo2{$x} += $y * $y;
	$histoc{$x}++;
}

foreach my $bin (sort {$a <=> $b} keys %histoc) {
	my $n = $histoc{$bin};
	my $mw = $histo{$bin} / $n;
	my $tmp = $histo2{$bin} / $n - $mw * $mw;
	my $sd = ($tmp > 0.) ? sqrt($tmp) / ($n - 1.) : 0.;
	print "$bin $mw $sd $n\n";
}

