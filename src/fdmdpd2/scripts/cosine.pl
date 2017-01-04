#!/usr/bin/perl -w
#
# This program takes an initial configuration as input and extends it by
# placing n*n identical copies of the same system in the output configuration.
# The molecules are randomized.
# 
# Syntax: extendsys.pl [n] < input.dat > output.dat
#
use strict;

my @L = ();
my $n = 0;
my $N = 0;
my $t = 0.;
my $count = 0;
my @COORD;

# read in original results
while(<>) {
	if(/^# L=([\d\.]+) ([\d\.]+) ([\d\.]+) n=(\d+) N=(\d+) t=([\d\.]+)/) {
		$L[0] = $1;
		$L[1] = $2;
		$L[2] = $3;
		$n = $4;
		$N = $5;
		$t = $6;

		print "# L=$L[0] $L[1] $L[2] n=$n N=$N t=0\n";
		next;
	}
	elsif(/^# v=/) {
		print;
		next;
	}
	else {
		my @line = split;
		my $i = int($count / $N);
		my $j = $count - $i * $N;
		my @tmp;
		
		foreach my $blub (@line) {
			push @tmp, $blub;
		}
		$tmp[0] -= $L[0] * (-.25 + .2 * cos($tmp[0] * 2. * 3.141592654 / $L[0]));
		my $new_line = join(' ' , @tmp) . "\n";
		$COORD[2*$i+0][$j] = $new_line;
	
		@tmp = ();	
		foreach my $blub (@line) {
			push @tmp, $blub;
		}
		$tmp[0] += $L[0] * (-.25 + .2 * cos($tmp[0] * 2. * 3.141592654 / $L[0]));
		$new_line = join(' ' , @tmp) . "\n";
		$COORD[2*$i+1][$j] = $new_line;

		$count++;
	}
}

fisher_yates_shuffle( \@COORD );

for(my $i=0;$i<2*$n;$i++) {
	for(my $j=0;$j<$N;$j++) {
		print $COORD[$i][$j];
	}
}

# fisher_yates_shuffle( \@array ) : generate a random permutation
# of @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

