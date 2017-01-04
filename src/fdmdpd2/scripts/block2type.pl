#!/usr/bin/perl -w
use strict;

my $block = -1;

while(<>) {
	if (/^# L=([^ ]+) ([^ ]+) ([^ ]+) t=([^ ]+) blocks=(\d+)/) {
		print;
                $block = -1;
	} elsif (/^# n=(\d+) N=(\d+) name=/) {
		print;
		$block++;
	} elsif (/^# /) {
		print;
	} else {
		chop;
		my @tmp = split / /;
		$tmp[6]=$block;
		print join(' ', @tmp) . "\n";
	}
}

