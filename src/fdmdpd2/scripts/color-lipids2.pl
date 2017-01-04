#!/usr/bin/perl -w
#
# Syntax: color-lipids.pl [N] [vtf-file]
# 
# N: uniform chain length
#
# Read in vtf file and produce many more vtf-files where the 
# orientational order paramter of the lipids is encoded as "charge".
#

use strict;
use Math::Trig;

sub syntax
{
	printf("Syntax: color-lipids.pl [N] [vtf-file]\n\n");
	exit(0);
}

my @n = (1., 0., 0.);

		
my $total_mean = 0.;
my $total_mean2 = 0.;
my $total_count = 0;


my $N = shift or syntax();
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
				print $line;
			}

			if ($line =~ /^b\w*/) {
				print $line;
			}

			if ($line =~ /^pbc/) {
				if ($coordcount > 0) {
					my @lambdas = ();
					calc_lambda($coordcount, $N, \@type, \@x, \@y, \@z, \@lambdas);
					my $n = $coordcount / $N;
					print STDERR "I have found $n molecules, $coordcount!\n";

					print "\ntimestep\n$pbcline";
					for (my $i = 0; $i < $n; $i++) {
						for (my $j = 0; $j < $N; $j++) {
							print $x[$i * $N + $j] . " ";
							print $y[$i * $N + $j] . " ";
							print $z[$i * $N + $j] . " ";
							print $lambdas[$i] . " 0 0\n";
						}
					}

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

		for(my $i = 0; $i < $N - 2; $i++) {
			if (!defined $$t[$i]) {
				print STDERR "Fehler: i=$i j=$j\n";
			} 
			elsif ($$t[$i] == 0 && $$t[$i + 1] == 0) {
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

