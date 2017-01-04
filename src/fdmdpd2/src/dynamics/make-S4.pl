#!/usr/bin/perl -w
use strict;

my @dp = ();
my $c = 0;
my $dq = 0.05;
my %histo1 = ();
my %histo2 = ();
my %histo3 = ();
my %histoc = ();

while(<>) {
	next if (/^\D/);
	# t q1 q2 |q| SS SD DD
	my ($t, $qy, $qz, $q, $ss, $sd, $dd) = split / /;
	my $bin = int($q / $dq);
	$histo1{$t}[$bin] += $ss;
	$histo2{$t}[$bin] += $sd;
	$histo3{$t}[$bin] += $dd;
	$histoc{$t}[$bin]++;
}

print "# t |q| SS SD DD\n"; 
foreach my $time (sort {$a <=> $b} (keys(%histoc))) {

	for (my $i = 0; $i < @{$histoc{$time}}; $i++) {
		my $q = $i * $dq;
		next if (!(defined $histoc{$time}[$i]));
		my $h1 = $histo1{$time}[$i] / $histoc{$time}[$i];
		my $h2 = $histo2{$time}[$i] / $histoc{$time}[$i];
		my $h3 = $histo3{$time}[$i] / $histoc{$time}[$i];
		my $q2t = $q * $q * $time;

		print "$time $q $h1 $h2 $h3\n";
	}
	print "\n";
}
