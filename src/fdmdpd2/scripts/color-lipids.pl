#!/usr/bin/perl -w
#
# Syntax: color-lipids.pl [N] [threshold] [vtf-file]
# 
# N: uniform chain length
# threshold: critical value for order parameter
#
# Read in vtf file and produce many more vtf-files where the lipids with an
# order-parameters above the threshold are especially marked.
#

use strict;
use Math::Trig;

sub syntax
{
	printf("Syntax: color-lipids.pl [N] [threshold] [vtf-file]\n\n");
	exit(0);
}

my @n = (1., 0., 0.);

		
my $total_mean = 0.;
my $total_mean2 = 0.;
my $total_count = 0;


my $N = shift or syntax();
my $threshold = shift or syntax();
my $outcount = 0;

while(my $filename = shift)
{
	if (-e $filename) {
		# printf("Opening $filename...\n");
		open(FH, "<$filename") or die "Cannot open $filename for reading: $!\n\n";

		my $time;
		my @L;
		my $coordcount = 0;
		my @atom;
		my @type;
		my @bondline;
		my $blc = 0;
		my @x;
		my @y;
		my @z;
		my $pbcline;

		while(my $line = <FH>) {

			if ($line =~ /^a\w* (\d+) (.*)$/) {
				my $idx = $1;
				$atom[$idx] = $line;
				$atom[$idx] =~ s/\n//g;
				$line =~ /n\w*\s(\w+)\s/;
				$type[$idx] = 0 if ($1 eq 'A');
				$type[$idx] = 1 if ($1 eq 'B');
			}

			if ($line =~ /^b\w*/) {
				$bondline[$blc++] = $line;
			}

			if ($line =~ /^pbc/) {
				if ($coordcount > 0) {
					my @lambdas = ();
					calc_lambda($coordcount, $N, \@type, \@x, \@y, \@z, \@lambdas);
					my $n = $coordcount / $N;
					print "blc=$blc\n";
					print "I have found $n molecules, $coordcount!\n";
					open(OUT, ">out.$outcount.vtf") or die "$!\n";

					for (my $i = 0; $i < $n; $i++) {
						for (my $j = 0; $j < $N; $j++) {
							print OUT $atom[$i * $N + $j];
							if ($lambdas[$i] > $threshold) {
								print OUT " s GEL\n";
							} else {
								print OUT " s FLU\n";
							}
						}
					}

					for (my $i = 0; $i < $blc; $i++) {
						print OUT $bondline[$i];
					}
					
					print OUT "\ntimestep\n$pbcline";
					for (my $i = 0; $i < $n; $i++) {
						for (my $j = 0; $j < $N; $j++) {
							print OUT $x[$i * $N + $j] . " ";
							print OUT $y[$i * $N + $j] . " ";
							print OUT $z[$i * $N + $j] . "\n";
						}
					}

					close(OUT);
					$outcount++;
				}
				$coordcount = 0;
				@x = ();
				@y = ();
				@z = ();
				$pbcline = $line;
			}

			if ($line =~ /^([\.\de\+\-]+) ([\.\de\+\-]+) ([\.\de\+\-]+)/) {
				$x[$coordcount] = $1;
				$y[$coordcount] = $2;
				$z[$coordcount] = $3;
				$coordcount++;
			}
		}
	}
}

sub calc_lambda
{
	my $nN = shift;
	my $N = shift;
	my $n = $nN / $N;
	my $t = shift;
	my $x = shift;
	my $y = shift;
	my $z = shift;
	my $lambdas = shift;
	
	for (my $j = 0; $j < $n; $j++) {
		my $mean = 0.;
		my $count = 0;

		for(my $i = 0; $i < $N - 1; $i++) {
			if ($$t[$i] == 0 && $$t[$i + 1] == 0) {
				my $dx = $$x[$j * $N + $i + 1] - $$x[$j * $N + $i];
				my $dy = $$y[$j * $N + $i + 1] - $$y[$j * $N + $i];
				my $dz = $$z[$j * $N + $i + 1] - $$z[$j * $N + $i];
				
				my $ct2 = $dx*$dx/($dx*$dx+$dy*$dy+$dz*$dz);
				my $sxx = 1.5 * $ct2 - .5;

				$mean += $sxx;
				$count++;
			}
		}
		$mean /= $count;
#		print "j=$j, mean=$mean, count=$count\n";
		push @{$lambdas}, $mean;
	}
}

