#!/usr/bin/perl -w
use strict;

my @data = ();
my @count = ();
my $cnt = 0;

while(my $fn = shift) {
	open(FH, "<$fn") or die "Couldn't open $fn for reading: $!\n";
	my $line = 0;
	$cnt++;

	while(<FH>) {
		next if(/^#/);
		chomp;
		my @tmp = split /\s/;

		for (my $i = 0; $i < @tmp; $i++) {
			$data[$line][$i] += $tmp[$i];
			$count[$line][$i]++;
		}
		$line++;
	}
	close(FH);
}
print STDERR "Averaging over $cnt graphs.\n";

for (my $i = 0; $i < @data; $i++) {
	for (my $j = 0; $j < @{$data[$i]}; $j++) {
		my $a = $data[$i][$j] / $count[$i][$j];
		print "$a";
		if ($j < @{$data[$i]} - 1) {
			print " ";
		} else {
			print "\n";
		}
#		if ($count[$i][$j] != $cnt) {
#			print STDERR "cnt=$cnt count=$count[$i][$j]\n\n";
#		}
	}
}
