#!/usr/bin/perl -w
#
# Perl script to make a histogram with error bars
#
# Syntax: ./mkhisto.pl [bw] {file}
#
use strict;

sub syntax
{
	print STDERR "Syntax: [bw] {file}\n\n";
	exit 0;
}

my $bw = shift or syntax;
my %histo = ();
my %histo2 = ();
my %histoc = ();

while(<>) {
	next if (/^#/);
	my @tmp = split / /;
	my $x = $tmp[0];
	my $y = $tmp[1];

	my $bin = int($x / $bw);

	$histo{$bin} += $y;
	$histo2{$bin} += $y * $y;
	$histoc{$bin}++;
}

foreach my $bin (sort {$a <=> $b} keys %histoc) {
	my $x = ($bin + .5) * $bw;
	my $n = $histoc{$bin};
	my $mw = $histo{$bin} / $n;
	my $tmp = $histo2{$bin} / $n - $mw * $mw;
	my $sd = ($tmp > 0.) ? sqrt($tmp) / ($n - 1.) : 0.;
	print "$x $mw $sd $n\n";
}

