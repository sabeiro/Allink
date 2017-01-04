#!/usr/bin/perl -w
#
# If an isolated peak at 0 shows up, then $delta (lateral binwidth) is too
# small.
#
use strict;
use POSIX;

my $delta = 2.0;
my $dx = 0.05;

my $dy;
my $dz;
my $xbins;
my $ybins;
my $zbins;

my @L;
my $blocks;
my $block;
my $k;
my $nN;
my $name;

my $files = 0;
my @storage;
my %grid;
my $minbin = 10000;
my $maxbin = 10000;
my @midpoint = ();
my @midpoint_c = ();


while(<>) {
	if (/^# L=([\d\.]+) ([\d\.]+) ([\d\.]+) t=([\d\.]+) blocks=(\d+)/) {
		print STDERR;
		@L = ($1, $2, $3);
		$blocks = $5;
		$files++;

		$xbins = int($L[0] / $dx);
		$ybins = int($L[1] / $delta) + 1;
		$zbins = int($L[2] / $delta) + 1;
		$dy = $L[1] / $ybins;
		$dz = $L[2] / $zbins;

		@storage = ();
		$block = 0;
		$nN = 0;
		$k = 0;
		@midpoint = ();
		@midpoint_c = ();
		
	} elsif (/^# n=(\d+) N=(\d+) name=(\w+)/) {
		$nN += $1 * $2;
		$name = $3;
		$block++;
	} elsif (/^# /) {

	} else {
		chop;
		my @tmp = split / /;

		for my $d (0 .. 2) {
			while ($tmp[$d] >= $L[$d]) { $tmp[$d] -= $L[$d]; }
			while ($tmp[$d] < 0.) { $tmp[$d] += $L[$d]; }
		}

		my $ybin = int($tmp[1] / $dy);
		my $zbin = int($tmp[2] / $dz);
		my $bin = $zbins * $ybin + $zbin;
		my $t = int($tmp[6]);

		$midpoint[$bin] += $tmp[0];
		$midpoint_c[$bin]++;

		$storage[$k]{'x'} = $tmp[0];
		$storage[$k]{'bin'} = $bin;
		$storage[$k]{'type'} = $t;
		$storage[$k]{'name'} = $name;
		$k++;

		if ($k == $nN && $block == $blocks) {
			for (my $i = 0; $i < $ybins * $zbins; $i++) {
				if (defined $midpoint_c[$i]) {
					$midpoint[$i] /= $midpoint_c[$i];
					if ($midpoint_c[$i] == 1) {
						warn("only 1 particle in bin $i\n");
					}
				}
			}
			for (my $i = 0; $i < $k; $i++) {
				my $bin = $storage[$i]{'bin'};
				my $x = $storage[$i]{'x'} - $midpoint[$bin];
				my $xbin = int(floor($x / $dx)) + 10000;
				my $t = $storage[$i]{'type'};
				$name = $storage[$i]{'name'};

				$minbin = $xbin if ($xbin < $minbin);
				$maxbin = $xbin if ($xbin > $maxbin);
				$grid{$name}[$t][$xbin] += 1./($L[1] * $L[2]);
			}
		}
	}
}

print "# ";
foreach $name (keys %grid) { 
	for (my $t = 0; $t < @{$grid{$name}}; $t++) {
		print "${name}_$t ";
	}
}
print "\n";

for (my $i = $minbin; $i < $maxbin; $i++) {
	my $x = ($i - 10000 + .5) * $dx;
	print "$x";
	foreach $name (keys %grid) { 
		for (my $t = 0; $t < @{$grid{$name}}; $t++) {
			my $r = 0;
			if (defined $grid{$name}[$t][$i]) {
				$r = $grid{$name}[$t][$i];
				$r /= $files * $dx;
			}
			print " $r";
		}
	}
	print "\n";
}

