#!/usr/bin/perl -w
use strict;

sub syntax
{
	print STDERR "Syntax: $0 [n2_min] [n2_max]\n\n";
	exit 1;
}


my $Kx = 128;
my $Ky = 128;

my $n2_min = shift or syntax();
my $n2_max = shift or syntax();

#my $n2_min = 4500; # inclusive 
#my $n2_max = 5000; #exclusive;

# split input files. each mode goes into its own file, this may not work
# everywhere, because many files will be opened at the same time
my @files;
for (my $i = 0; $i < $Kx; $i++) {
	my $nx = ($i < $Kx / 2) ? $i : $i - $Kx;
	for (my $j = 0; $j < $Ky / 2 + 1; $j++) {
		my $ny = $j;

		my $n2 = $nx * $nx + $ny * $ny;
		if ($n2 >= $n2_min && $n2 < $n2_max) {
			open($files[$i][$j], ">>u_${nx}_${ny}") or die "Opening: $nx $ny: $!";
		}
	}
}

while(<>) { # consume header line 
	for (my $i = $Kx / 2; $i < $Kx; $i++) {
		my $nx = ($i < $Kx / 2) ? $i : $i - $Kx;
		for (my $j = 0; $j < $Ky / 2 + 1; $j++) {
			my $ny = $j;
			my $line = <>;

			my $n2 = $nx * $nx + $ny * $ny;
			if ($n2 >= $n2_min && $n2 < $n2_max) {
				my $fh = $files[$i][$j];
				print $fh "$line";
			}
		}
#		<>;
	}
	for (my $i = 0; $i < $Kx / 2; $i++) {
		my $nx = $i; #($i < $Kx / 2) ? $i : $i - $Kx;
		for (my $j = 0; $j < $Ky / 2 + 1; $j++) {
#			next if ($i == 0 && $j == 0);
			my $ny = $j;
			my $line = <>;

			my $n2 = $nx * $nx + $ny * $ny;
			if ($n2 >= $n2_min && $n2 < $n2_max) {
				my $fh = $files[$i][$j];
				print $fh "$line";
			}
		}
#		<> if ($i < $Kx / 2 - 1);
	}
}

for (my $i = 0; $i < $Kx; $i++) {
	my $nx = ($i < $Kx / 2) ? $i : $i - $Kx;
	for (my $j = 0; $j < $Ky / 2 + 1; $j++) {
		my $ny = $j;
		my $n2 = $nx * $nx + $ny * $ny;
		if ($n2 >= $n2_min && $n2 < $n2_max) {
			close($files[$i][$j]) or die "Closing: $i $j: $!";
		}
	}
}
