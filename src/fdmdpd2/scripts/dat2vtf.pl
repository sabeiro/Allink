#!/usr/bin/perl -w
# 
# dat2vtf.pl - Read from STDIN and write extended VTF file to STDOUT
#
# "extended" means, that the coordinates and also the velocities are printed
# The original vtf plugin from VMD ignores the velocities, however the modified
# VTF plugin understands them.
#
use strict;
use POSIX qw/floor/;

# particle properties, add more if needed.
my @names = qw/A B C/;
my @radii = qw/0.8 0.8 0.8/;

my $print_header = 1;
my @L;
my $time;
my $blocks;
my $block;
my $nN;
my $c;
my @n;
my @N;
my @resname;
my @coord;
my @type;

while(<>) {
	chomp;
	if (/^# L=([^ ]+) ([^ ]+) ([^ ]+) t=([^ ]+) blocks=(\d+)/) {
		print STDERR "$_\n";
		$L[0] = $1;
		$L[1] = $2;
		$L[2] = $3;
		$time = $4;
		$blocks = $5;
		$nN = 0;
		$block = 0;
		$c = 0;
	}
	elsif (/^# n=(\d+) N=(\d+) name=(\w+)/) {
		$nN += $1 * $2;
		$n[$block] = $1;
		$N[$block] = $2;
		$resname[$block] = $3;
		$block++;
	}
	elsif (/^#/) { 

	}
	else {
		my @tmp = split / /;
		for (my $i = 0; $i < 6; $i++) {
			$coord[$c][$i] = shift @tmp;
		}
		$type[$c++] = shift @tmp;
	}

	if ($c == $nN && $block == $blocks) {
		# print header if necessary
		if ($print_header) {
			my $chain = 0;
			my $idx = 0;
			for (my $i = 0; $i < $blocks; $i++) {
				for (my $j = 0; $j < $n[$i]; $j++, $chain++) {
					my $idx0 = $idx;
					for (my $k = 0; $k < $N[$i]; $k++, $idx++) {
						my $name = $names[$type[$idx]];
						my $radius = $radii[$type[$idx]];
						printf("a %llu r %lg n %s resid %d res %s\n",
							$idx, $radius, $name,
							$chain, $resname[$i]);
					}
					printf("b %d::%d\n", $idx0, $idx - 1);
				}
			}
			$print_header = 0;
		}

		# print coordinates **and** velocities(!)
		printf("\ntimestep\npbc %lg %lg %lg\n", $L[0], $L[1], $L[2]);
		my $idx = 0;
		for (my $i = 0; $i < $blocks; $i++) {
			for (my $j = 0; $j < $n[$i]; $j++) {
				my @xv_com;
				my $idx0 = $idx;

				# calculate center-of-mass
				for (my $k = 0; $k < $N[$i]; $k++, $idx++) {
					for (my $d = 0; $d < 3; $d++) {
						$xv_com[$d] += $coord[$idx][$d];
					}
				}
				# shift center-of-mass, so that it is in the primary box
				for (my $d = 0; $d < 3; $d++) {
					$xv_com[$d] /= $N[$i];
					$xv_com[$d] = $L[$d] * floor($xv_com[$d]/ $L[$d]);
				}
				for (my $k = 0, $idx = $idx0; $k < $N[$i]; $k++, $idx++) {
					for (my $d = 0; $d < 3; $d++) {
						$coord[$idx][$d] -= $xv_com[$d];
					}
					printf("%lg %lg %lg %lg %lg %lg\n",
						$coord[$idx][0], $coord[$idx][1],
						$coord[$idx][2], $coord[$idx][3],
						$coord[$idx][4], $coord[$idx][5]);
				}
			}
		}
	}
}

