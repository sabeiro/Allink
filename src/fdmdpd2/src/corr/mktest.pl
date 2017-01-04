#!/usr/bin/perl -w
use strict;

for (my $i = 0; $i < 256; $i++) {
	my $a = exp(-.01*$i);
	print "$i $a\n";
}
