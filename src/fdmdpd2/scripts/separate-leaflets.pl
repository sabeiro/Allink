#!/usr/bin/perl -w
use strict;

my $n;
my $N;

my @storage;
my $xa = 0.;
my $xb = 0.;
my $xan = 0;
my $xbn = 0;

my @storage_u;
my @storage_l;
my $un = 0;
my $ln = 0;

while(<>) {
	chop;
	if (/^#/) {
		s/^# //;
		if (/L=/) {
			s/blocks=\d+/blocks=2/;
			print "# $_\n";
		} elsif (/n=/) {
			s/\w+=//g;
			my @tmp = split / /;
			$n = shift @tmp;
			$N = shift @tmp;
		} else {
			print "# $_\n";
		}
	} else {
		push @storage, $_;
		my @tmp = split / /;
		if ($tmp[6] == 0) {
			$xa += $tmp[0];
			$xan++;
		} else {
			$xb += $tmp[0];
			$xbn++;
		}

		if ($xan + $xbn == $N) {
			if ($xb / $xbn > $xa / $xan) {
				# upper leaflet
				for (my $i = 0; $i < $N; $i++) {
					push @storage_u, $storage[$i] . "\n";
				}
				$un++;
			} else {
				# lower leaflet
				for (my $i = 0; $i < $N; $i++) {
					push @storage_l, $storage[$i] . "\n";
				}
				$ln++;
			}
			$xan = 0.;
			$xa = 0.;
			$xbn = 0.;
			$xb = 0.;
			@storage = ();
		}
	
		if ($un + $ln == $n) {
			print "# n=$un N=$N name=UPPER\n";
			for (my $i = 0; $i < $un * $N; $i++) {
				print $storage_u[$i];
			}
		
			print "# n=$ln N=$N name=LOWER\n";
			for (my $i = 0; $i < $ln * $N; $i++) {
				print $storage_l[$i];
			}
			@storage_l = ();
			@storage_u = ();
			$ln = 0;
			$un = 0;
		}
	}
}

