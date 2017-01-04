#!/usr/bin/perl -w
use strict;

my @dp = ();
my $c = 0;
my $dq = 0.05;
my %histo = ();
my %histoc = ();

while(<>) {
	next if (/^\D/);
	# t q_y q_z q Re(Fd) Im(Fd) Re(Fs) Im(Fs) Fsg(q,t)
	my ($t, $qy, $qz, $q, $re, $im, $G, $Y) = split / /;
	my $bin = int($q / $dq);
	if ($Y > 0.) {
		$histo{$t}[$bin] += ($re - $Y) / $Y;
		$histoc{$t}[$bin]++;
	}
}

print "# t |q| Y\n"; 
foreach my $time (sort {$a <=> $b} (keys(%histo))) {

	for (my $i = 0; $i < @{$histo{$time}}; $i++) {
		my $q = $i * $dq;
		next if (!(defined $histo{$time}[$i]));
		my $h = $histo{$time}[$i] / $histoc{$time}[$i];
		my $q2t = $q * $q * $time;

		print "$time $q $h\n";
	}
	print "\n";
}
