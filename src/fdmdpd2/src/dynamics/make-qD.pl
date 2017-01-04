#!/usr/bin/perl -w
use strict;

my @dp = ();
my $c = 0;
my $dq = 0.05;
my @histo = ();
my @histoc = ();

while(<>) {
	next if (/^\D/);
	# t q_y q_z q Re(Fd) Im(Fd) Re(Fs) Im(Fs)
	my ($t, $qy, $qz, $q, $rf, $if, $rs, $is) = split / /;
	$dp[$c]{'t'} = $t;
	$dp[$c]{'k'} = $q;
	$dp[$c]{'F'} = $rs;
	$c++;

	my $bin = int($q / $dq);
	$histo[$bin]{$t} += $rs;
	$histoc[$bin]{$t}++;
}

for (my $i = 0; $i < @histo; $i++) {
	my $q = $i * $dq;
	next if (!(defined $histo[$i]));
#	print STDERR "$q \n";
	my %hash = %{$histo[$i]};

	open(FILE, ">data.tmp") or die "Cannot open data.tmp for writing: $!\n";

	foreach my $time (sort {$a <=> $b} (keys(%{$histo[$i]}))) {
		my $h = $histo[$i]{$time} / $histoc[$i]{$time};
		my $q2t = $q * $q * $time;
#		print FILE "$q $time $q2t $h\n";
		print FILE "$q2t $h\n";
#		print "$q2t $h\n" if ($h > 0.);
	}

	close(FILE);

	open(GNUPLOT, "|gnuplot") or die "Cannot open gnuplot for pipe: $!\n";
	print GNUPLOT "D=1e-5\n";
	print GNUPLOT "fit exp(-D*x) \"data.tmp\" u 1:2 via D\n";
	close(GNUPLOT);

	open(FIT, "<fit.log") or die "Cannot open fit.log: $!\n";

	my $D = -1.;
	while(<FIT>) {
		if(/^D\s+= ([^ ]+)\s+\+\/\-/) {
			$D = $1;
		}
	}

	close(FIT);
	unlink('fit.log');

	print "$q $D\n";
}
