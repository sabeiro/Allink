#!/usr/bin/perl -w
# velocity profile calculation for RNEMD method

$dg = 1;	# dimension of velocity gradient
$dv = 2;	# dimension of velocity direction

$bw = 0.2;

$files = 0;

while(<>) {
	if (/^# L=([^ ]+) ([^ ]+) ([^ ]+)/) {
		my @len = ($1, $2, $3);
		for ($i = 0, $A = 1.; $i < 3; $i++) {
			if ($i == $dg) {
				$L = $len[$i];
			} else {
				$A *= $len[$i];
			}
		}
		$L = $len[$dg];
		$files++;
	} elsif (/^#/) {

	} else {
		my @tmp = split / /;

		# position
		while ($tmp[$dg] >= $L) { $tmp[$dg] -= $L; }
		while ($tmp[$dg] < 0.) { $tmp[$dg] += $L; }
		$bin = int($tmp[$dg] / $bw);
	
		# velocity
		$v = $tmp[3 + $dv];
		$histo_v[$bin] += $v;
		$histo_v2[$bin] += $v ** 2;
	
		# temperature
		$T = ($tmp[3] ** 2 + $tmp[4] ** 2 + $tmp[5] ** 2) / 3.;
		$histo_T[$bin] += $T;
		$histo_T2[$bin] += $T ** 2;
	
		# density / normalization
		$histoc[$bin]++;
	}
}

print "# y v_z D(v_z) T D(T) rho D(rho) n\n";
for (my $i = 0; $i < @histo_v; $i++) {
	$n = $histoc[$i];
	$x = ($i + .5) * $bw;
	print "$x ";

	# velocity
	$mw = $histo_v[$i] / $n;
	$sd = sqrt(($histo_v2[$i] / $n - $mw * $mw) / ($n - 1.));
	print "$mw $sd ";

	# temperature
	$T = $histo_T[$i] / $n;
	$sd = sqrt(($histo_T2[$i] / $n - $T * $T) / ($n - 1.));
	print "$T $sd ";

	# density
	$rho = $n / $files / $A;
	print "$rho $n\n";
}
