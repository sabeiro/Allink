#!/usr/bin/perl -w
use strict;

my $n = shift;
my $start = shift;
syntax() if (!defined($n) || $n == 0);
syntax() if (!defined($start)  || $start == 0);

my $cnt = 0;
my @data = ();
my $i0 = 0;

while(my $line = <>) {
	next if ($line =~ /^# /);
	chomp $line;

	for (my $i = $i0; $i < 1 + int($cnt / $start); $i++) {
		if ($cnt - $i * $start < $n) {
			push @{$data[$i]}, $line;
		}
	}

	if (@{$data[$i0]} == $n) {
		open(FH, sprintf(">out.%05d", $i0)) or die "Couldn't open out.$i0 for writing: $!";
		for (my $j = 0; $j < @{$data[$i0]}; $j++) {
			my $a = $data[$i0][$j];
			print FH "$a\n";
		}
		@{$data[$i0]} = ();
		$i0++;
	}

	$cnt++;
}

for (my $i = 0; $i < @data; $i++) {
	next if (@{$data[$i]} < $n);

	open(FH, sprintf(">out.%05d", $i)) or die "Couldn't open out.$i for writing: $!";
	for (my $j = 0; $j < @{$data[$i]}; $j++) {
		my $a = $data[$i][$j];
		print FH "$a\n";
	}
	close(FH);
}

sub syntax
{
	die "Syntax: $0 [window size] [step size]\n\n";
	
}


