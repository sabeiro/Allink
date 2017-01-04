#!/usr/bin/perl -w
use strict;

my $dr = 0.01;
my @L;
my $n;
my $N;
my $name;

my @x;
my %g;
my $files = 0;

while(<>) {
	if (/^# L=([^ ]+) ([^ ]+) ([^ ]+)/) {
		$L[0] = $1;
		$L[1] = $2;
		$L[2] = $3;

		$files++;
		print STDERR;
	} elsif (/^# n=(\d+) N=(\d+) name=(\w+)/) {
		$n = $1;
		$N = $2;
		$name = $3;

		@x = ();
	} elsif (/^#/) {

	} else {
		my @tmp = split / /;
		push @x, { 'x' => $tmp[0], 'y' => $tmp[1], 'z' => $tmp[2],
			'vx' => $tmp[3], 'vy' => $tmp[4], 'vz' => $tmp[5],
			't' => $tmp[6] };

		if (@x == $n * $N) {
			for my $i (0 .. $#x) {
#				print "$i / $#x\n";
				for my $j ($i+1 .. $#x) {
					my $dx = $x[$i]{'x'} - $x[$j]{'x'};
					my $dy = $x[$i]{'y'} - $x[$j]{'y'};
					my $dz = $x[$i]{'z'} - $x[$j]{'z'};

#					while ($dx > .5 * $L[0]) { $dx -= $L[0]; }
#					while ($dx < -.5 * $L[0]) { $dx += $L[0]; }
					while ($dy > .5 * $L[1]) { $dy -= $L[1]; }
					while ($dy < -.5 * $L[1]) { $dy += $L[1]; }
					while ($dz > .5 * $L[2]) { $dz -= $L[2]; }
					while ($dz < -.5 * $L[2]) { $dz += $L[2]; }
					
					my $r = sqrt($dy**2 + $dz**2);
					my $bin = int($r / $dr);
					$g{$name}[$bin]++;
				}
			}
		}
	}
}

print "# r ";
my $max = 0;
foreach $name (keys %g) {
	print "$name ";
	$max = $#{$g{$name}} if ($#{$g{$name}} > $max);
}
print "\n";

my $pref = $L[1] * $L[2] / ($files * $n * $n * $N * $N);
for my $i (0 .. $max) {
	my $l = $i * $dr;
#	my $r = ($l**3 + 1.5*$l*$l*$dr + 1.5*$l*$dr*$dr + .5*$dr*$dr*$dr) ** (1./3.);
	my $r = sqrt($l * $l + $l * $dr + .5 * $dr * $dr);
	my $dv = 2. * 3.141592654 * $r * $dr;
	print "", ($i + .5) * $dr, " ";
	foreach $name (keys %g) {
		my $g = 0.;
		if (defined $g{$name}[$i]) {
			$g = 2. * $g{$name}[$i] / $dv *  $pref;
		}	
		print "$g ";
	}
	print "\n";	
}

