#!/usr/bin/perl -w
use strict;

my $na = shift;
my $nb = shift;
my @arch;

for(my $i = 0; $i < $na; $i++) {
	$arch[$i] = 0;
}
for(my $i = $na; $i < $na+$nb; $i++) {
	$arch[$i] = 1;
}

my $count = 0;

while(<>) {
	if (!(/^#/)) {
		my $bead = $arch[$count++ % ($#arch+1)];
		s/(\d+)$/$bead/;
	}
	print;
}
