#!/usr/bin/perl -w
use strict;

my $Kx = 16;
my $Ky = 16;

# split input files. each mode goes into its own file, this may not work
# everywhere, because many files will be opened at the same time
my @files;
for (my $i = 0; $i < $Kx; $i++) {
	my $nx = ($i < $Kx / 2) ? $i : $i - $Kx;
	for (my $j = 0; $j < $Ky / 2 + 1; $j++) {
		my $ny = $j;

		open($files[$i][$j], ">u_${nx}_${ny}") or die "Opening: $nx $ny: $!";
	}
}

while(<>) { # consume header line 
	for (my $i = 0; $i < $Kx; $i++) {
		for (my $j = 0; $j < $Ky / 2 + 1; $j++) {
			my $line = <>;
			my $fh = $files[$i][$j];
			print $fh "$line";
		}
	}
}

for (my $i = 0; $i < $Kx; $i++) {
	for (my $j = 0; $j < $Ky / 2 + 1; $j++) {
		close($files[$i][$j]) or die "Closing: $i $j: $!";
	}
}
