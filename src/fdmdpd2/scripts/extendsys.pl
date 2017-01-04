#!/usr/bin/perl -w
#
# This program takes an initial configuration as input and extends it by
# placing identical copies of the same system in the output configuration.
# The molecules are randomized.
# 
# Syntax: extendsys.pl [nx] [ny] [nz] < input.dat > output.dat
#
use strict;

sub syntax
{
	print "Syntax: $0 [nx] [ny] [nz]\n\n";
	exit 0;
}

my $factor_nx = int(shift) or syntax();
my $factor_ny = int(shift) or syntax();
my $factor_nz = int(shift) or syntax();
my @L = ();
my $n = 0;
my $N = 0;
my $blocks = 0;
my $block = 0;
my $count = 0;
my @COORD = ( );

# read in original results
while(<>) {
	chomp;
	if (/^# L=([^ ]+) ([^ ]+) ([^ ]+) t=([^ ]+) blocks=(\d+)/) {
		$L[0] = $1;
		$L[1] = $2;
		$L[2] = $3;
		my $time = $4;
		$blocks = $5;

		my $new_L1 = $L[0] * $factor_nx;
		my $new_L2 = $L[1] * $factor_ny;
		my $new_L3 = $L[2] * $factor_nz;
		print "# L=$new_L1 $new_L2 $new_L3 t=$time blocks=$blocks\n";

	} elsif (/^# n=([^ ]+) N=([^ ]+) name=([^ ]+)/) {
		$n = $1;
		$N = $2;
		my $name = $3;

		my $new_n = $n * $factor_nx * $factor_ny * $factor_nz;
		print "# n=$new_n N=$N name=$name\n";
	} elsif (/^#/) {
		print;
		print "\n";

	} else {
		my @line = split;
		my $i = int($count / $N);
		my $j = $count - $i * $N;

		for (my $x=0;$x<$factor_nx;$x++) {
			for(my $y=0;$y<$factor_ny;$y++) {
				for(my $z=0;$z<$factor_nz;$z++) {
					my @tmp;
					foreach my $blub (@line) {
						push @tmp, $blub;
					}
					$tmp[0] += $L[0] * $x;
					$tmp[1] += $L[1] * $y;
					$tmp[2] += $L[2] * $z;
					my $new_line = join(' ' , @tmp) . "\n";
					$COORD[(($i * $factor_nx + $x) * $factor_ny + $y) * $factor_nz + $z][$j] = $new_line;
				}
			}
		}

		if (++$count == $n * $N) {
			fisher_yates_shuffle( \@COORD );

			for (my $i = 0; $i < $factor_nx * $factor_ny * $factor_nz * $n; $i++) {
				for(my $j = 0; $j < $N; $j++) {
					print $COORD[$i][$j];
				}
			}
			@COORD = ( );
			$count = 0;
		}
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

