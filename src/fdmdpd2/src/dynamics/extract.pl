#!/usr/bin/perl -w
#
# Syntax: extract.pl [resume.dat] [resume-isf.dat] [-g|-c|-s|-d]
#
# -g: extract static collective scattering function
# -c: extract collective intermediate scattering function
# -s: extract self intermediate scattering function
# -d: extract mean-square displacement, non-gaussianity
# -q: extract dynamic overlap parameter and its susceptibility
# -S4: extract S_4^{ol}(q,t) SS / SD / DD structure function
# -v: extract velocity autocorrelation function
# -F: extract force autocorrelation function
# -u: extract undulation autocorrelation function
# -U: extract undulation power spectrum
# -jt: extract autocorrelation of transversal current
# -jl: extract autocorrelation of longitudinal current
#
# The system file is only needed to get the box lengths.
#
#
use strict;
use constant PI => 4 * atan2(1, 1);
use List::Util qw[min max];


my $sysfilename = shift;
my $isffilename = shift;
my $option = shift;
my @L;
my $line;
my $mode;
my @G;
my @S0;
my @r2;
my @r4;
my @fp_Q;
my @S4;
my @va;
my @facf_g;
my @facf_h;
my @Cls0;
my @Cld0;
my @Cts0;
my @Ctd0;

my $isf_enabled = 0;
my $fp_enabled = 0;
my $vacf_enabled = 0;
my $facf_enabled = 0;
my $oldxv_enabled = 0;
my $uacf_enabled = 0;
my $jacf_enabled = 0;

if ($option eq '-g') {
	$mode = 1;
} elsif ($option eq '-c') {
	$mode = 2;
} elsif ($option eq '-s') {
	$mode = 3;
} elsif ($option eq '-d') {
	$mode = 4;
} elsif ($option eq '-q') {
	$mode = 5;
} elsif ($option eq '-S4') {
	$mode = 6;
} elsif ($option eq '-v') {
	$mode = 7;
} elsif ($option eq '-F') {
	$mode = 8;
} elsif ($option eq '-u') {
	$mode = 9;
} elsif ($option eq '-U') {
	$mode = 10;
} elsif ($option eq '-j') {
	$mode = 11;
} elsif ($option eq '-J') {
	$mode = 12;
} else {
	die "Option $option not understood. Goodbye.\n\n";
}


open(SYS, "<$sysfilename") or die "Cannot open $sysfilename: $!\n";
$line = <SYS>;
$line =~ /L=([^ ]+) ([^ ]+) ([^ ]+)/;
$L[0] = $1;
$L[1] = $2;
$L[2] = $3;
close(SYS);

# blocks=10 elements=19 count=2501 nN=4680 dt=0.02 obs=fs
# K=128 128
# a=0.3

open(ISF, "<$isffilename") or die "Cannot open $isffilename: $!\n";
$line = <ISF>;
chomp $line;
$line =~ /^# blocks=([^ ]+) elements=([^ ]+) count=([^ ]+) nN=([^ ]+) dt=([^ ]+) obs=(\w+)/;
my $max_blocks = $1;
my $max_elements = $2;
my $count = $3;
my $nN = $4;
my $dt = $5;
my $obs = $6;
print STDERR "debug: blocks=$max_blocks elements=$max_elements count=$count\n";
print STDERR "debug: nN=$nN dt=$dt obs=$obs\n";

my $K0 = -1;
my $K1 = -1;

# isf calculation
if ($obs =~ /[iuj]/) {
	$line = <ISF>;
	chomp $line;
	$line =~ /^# K=([^ ]+) ([^ ]+)/;
	$K0 = $1;
	$K1 = $2;
	print STDERR "debug: K=$K0 $K1\n";

	if ($obs =~ /i/) {
		$isf_enabled = 1;
		$oldxv_enabled = 1;
	}
	if ($obs =~ /u/) {
		$uacf_enabled = 1;
	}
	if ($obs =~ /j/) {
		$jacf_enabled = 1;
	}
}

# four-point / dynamic overlap "Q4"
if ($obs =~ /f/) {
	$line = <ISF>;
	chomp $line;
	$line =~ /^# a=([^ ]+)/;
	print STDERR "debug: a=$1\n";
	$fp_enabled = 1;
	$oldxv_enabled = 1;
}

# velocity autocorrelation
if ($obs =~ /v/) {
	$vacf_enabled = 1;
	$oldxv_enabled = 1;
}

# Force autocorrelation function
my $facf_dr = 0;
my $facf_rbins = 0;
if ($obs =~/F/) {
	$line = <ISF>;
	$line = <ISF>;
	chomp $line;
	$line =~ /^# facf_dr=([^ ]+) rbins=([^ ]+)/;
	$facf_dr = $1;
	$facf_rbins = $2;
	print STDERR "debug: rc=$1 facf_rbins=$2\n";
	$facf_enabled = 1;
	$oldxv_enabled = 1;
}

if (($mode == 1 || $mode == 2 || $mode == 3 || $mode == 4) && !($isf_enabled)) {
	print STDERR "mode=$mode, obs=$obs. No isf data included in file!\n";
	exit 1;
}

if (($mode == 5 || $mode == 6) && !($fp_enabled)) {
	print STDERR "mode=$mode, obs=$obs. No four-point data in file!\n";
	exit 1;
}

if (($mode == 7) && !($vacf_enabled)) {
	print STDERR "No VACF data found in this file!\n";
	exit 1;
}

if (($mode == 8) && !($facf_enabled)) {
	print STDERR "No FACF data found in this file!\n";
	exit 1;
}

if (($mode == 9) && !($uacf_enabled)) {
	print STDERR "No UACF data found in this file!\n";
	exit 1;
}

if (($mode == 11) && !($jacf_enabled)) {
	print STDERR "No JACF data found in this file!\n";
	exit 1;
}

my @blocklen;
my @c;

for (my $i = 0; $i < $max_blocks; $i++) {
	$blocklen[$i] = 0;
}
for (my $i = 0; $i < $max_blocks * $max_elements; $i++) {
	$c[$i] = 0;
}

my $blocks;
for (my $k = 1; $k < $count; $k++) {
	$blocks = 1;

	for (my $i = int($k / $max_elements); int($i) > 0; $i = int($i / $max_elements)) {
		$blocks++;
	}

	for (my $i = 0; $i < $blocks; $i++) {
		if ($k % int(($max_elements + 1) ** ($i) + .5) == 0) {
			#print "$k $i\n" if ($i > 0);
			$blocklen[$i]++;
			my $len = min($blocklen[$i], $max_elements);
			for (my $j = 0; $j < $len; $j++) {
				$c[$i * $max_elements + $j]++;
			}
		}
	}
}

if ($isf_enabled) {
	# read static collective structure factor
	for (my	$i = 0; $i < $K0 * ($K1 / 2 + 1); $i++) {
		$line = <ISF>;
		$line =~ /([\d\.\-\+e]+)/;
		$G[$i] = $1;
		$G[$i] /= $count * $nN;
	}
}

if ($uacf_enabled) {
	# read static undulation spectrum
	for (my	$i = 0; $i < $K0 * ($K1 / 2 + 1); $i++) {
		$line = <ISF>;
		$line =~ /([\d\.\-\+e]+)/;
		$S0[$i] = $1;
		$S0[$i] /= $count;
	}
}

if ($jacf_enabled) {
	# read static current correlation tensor
	my @q;

	for (my $n = 0; $n < $K0; $n++) {
		for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
			$line = <ISF>;
			$line =~ s/[ \n]+$//;
			my @tmp = split / /, $line;
			next if ($n == 0 && $m == 0);
			
			$q[0] = ($n < $K0 / 2) ? $n / $L[1] : ($n - $K0) / $L[1];
			$q[1] = $m / $L[2];
			my $inorm = 1. / sqrt($q[0] * $q[0] + $q[1] * $q[1]);
			$q[0] *= $inorm;
			$q[1] *= $inorm;

			my $k = $n * ($K1 / 2 + 1) + $m;
			$Cls0[$k] = 0.;
			$Cld0[$k] = 0.;
			$Cts0[$k] = 0.;
			$Ctd0[$k] = 0.;

			for my $ti (0 .. 1) {
				for my $tj (0 .. 1) {
					my $sr = shift @tmp;
					my $dr = shift @tmp;
					$sr /= $count;
					$dr /= $count;

					$Cls0[$k] += $q[$ti] * $q[$tj] * $sr;
					$Cld0[$k] += $q[$ti] * $q[$tj] * $dr;
					$Cts0[$k] -= $q[$ti] * $q[$tj] * $sr;
					$Ctd0[$k] -= $q[$ti] * $q[$tj] * $dr;

					if ($ti == $tj) {
						$Cts0[$k] += $sr;
						$Ctd0[$k] += $dr;
					}
				}
			}

			my $q1 = 2. * PI * (($n < $K0 / 2) ? $n / $L[1] : ($n - $K0) / $L[1]);
			my $q2 = 2. * PI * $m / $L[2];

			$Cls0[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
			$Cld0[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
			$Cts0[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
			$Ctd0[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
		}
	}
}

for (my $i = 0; $i < $blocks; $i++) {
	my $len = min($blocklen[$i], $max_elements);
	for (my $j = 0; $j < $len; $j++) {
		my @Sr = ();
		my @Si = ();
		my $time = ($j + 1) * $dt * (($max_elements + 1) ** $i);

		$line = <ISF>; # comment
		if (!($line =~ /block/)) {
			die "Achtung: header='$line'\n\n";
		}	
			
		if ($oldxv_enabled) {
			for (my $k = 0; $k < $nN; $k++) {
				<ISF>; # bead coordinates - ignored
			}
		}

		if ($isf_enabled) {
			$line = <ISF>; # r2 r4
			$line =~ /([\d\.\-\+e]+) ([\d\.\-\+e]+)/;
			$r2[$i * $max_elements + $j] = $1;
			$r4[$i * $max_elements + $j] = $2;

			for (my $k = 0; $k < $K0 * ($K1 / 2 + 1); $k++) {
				$line = <ISF>;
				if ($c[$i * $max_elements + $j] > 0 && ($mode == 2 || $mode == 3)) {
					$line =~ s/\n//;
					my @tmp = split(/ /, $line);
					if ($mode == 2) { # collective isf
						$Sr[$k] = $tmp[4];
						$Si[$k] = $tmp[5];
					} else { # single isf
						$Sr[$k] = $tmp[2];
						$Si[$k] = $tmp[3];
					}
					$Sr[$k] /= $nN;
					$Sr[$k] /= $c[$i * $max_elements + $j];
					$Si[$k] /= $nN;
					$Si[$k] /= $c[$i * $max_elements + $j];
				}
			}
	
			if ($c[$i * $max_elements + $j] > 0) {
				my $yr2 = $r2[$i * $max_elements + $j]; # standard msd
			       	$yr2 /= $nN * $c[$i * $max_elements + $j];
	
				if ($mode == 2 || $mode == 3) {
					print "# t q1 q2 |q| Re Im G Fsg\n";
		
					for (my $n = $K0 / 2; $n < $K0; $n++) {
						my $q1 = 2. * PI / $L[1] * ($n - $K0);
						for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
							my $q2 = 2. * PI / $L[2] * $m;
							my $q = sqrt($q1 * $q1 + $q2 * $q2);
							my $re = $Sr[$n * ($K1 / 2 + 1) + $m];
							my $im = $Si[$n * ($K1 / 2 + 1) + $m];
							my $g = $G[$n * ($K1 / 2 + 1) + $m];
							my $y = exp(-.25 * $q * $q * $yr2);
	
							print "$time $q1 $q2 $q $re $im $g $y\n";
						}
						print "\n";
					}
		
					for (my $n = 0; $n < $K0 / 2; $n++) {
						my $q1 = 2. * PI / $L[1] * $n;
						for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
							my $q2 = 2. * PI / $L[2] * $m;
							my $q = sqrt($q1 * $q1 + $q2 * $q2);
							my $re = $Sr[$n * ($K1 / 2 + 1) + $m];
							my $im = $Si[$n * ($K1 / 2 + 1) + $m];
							my $g = $G[$n * ($K1 / 2 + 1) + $m];
							my $y = exp(-.25 * $q * $q * $yr2);
	
							next if (($n == 0) && ($m == 0));
							print "$time $q1 $q2 $q $re $im $g $y\n";
						}
						print "\n";
					}
				}
			}
		}

		if ($fp_enabled) {
			$line = <ISF>; # <Q_S> <Q_D> <Q_S*Q_S> <Q_S*Q_D> <Q_D*Q_D>
			chomp $line;	
			my @tmp = split / /, $line;
			if ($#tmp < 4) {
				printf("Strange: ##$line##\n");
			}
			for (my $k = 0; $k < @tmp; $k++) {
				$fp_Q[$i * $max_elements + $j][$k] = $tmp[$k];
			}
	
			# read in S4 dynamic structure function SS / SD / DD
			for (my $k = 0; $k < $K0 * ($K1 / 2 + 1); $k++) {
				$line = <ISF>;
				if ($c[$i * $max_elements + $j] > 0 && ($mode == 6)) {
					$line =~ s/\n//;
					my @tmp = split(/ /, $line);
	
					$S4[3 * $k + 0] = $tmp[0];
					$S4[3 * $k + 0] /= $nN * $c[$i * $max_elements + $j];
					$S4[3 * $k + 1] = $tmp[1];
					$S4[3 * $k + 1] /= $nN * $c[$i * $max_elements + $j];
					$S4[3 * $k + 2] = $tmp[2];
					$S4[3 * $k + 2] /= $nN * $c[$i * $max_elements + $j];
				}
			}
	
			if ($c[$i * $max_elements + $j] > 0 && $mode == 6) {
				print "# t q1 q2 |q| SS SD DD\n";
		
				for (my $n = $K0 / 2; $n < $K0; $n++) {
					my $q1 = 2. * PI / $L[1] * ($n - $K0);
					for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
						my $q2 = 2. * PI / $L[2] * $m;
						my $q = sqrt($q1 * $q1 + $q2 * $q2);
						my $SS = $S4[3 * ($n * ($K1 / 2 + 1) + $m) + 0];
						my $SD = $S4[3 * ($n * ($K1 / 2 + 1) + $m) + 1];
						my $DD = $S4[3 * ($n * ($K1 / 2 + 1) + $m) + 2];
	
						print "$time $q1 $q2 $q $SS $SD $DD\n";
					}
					print "\n";
				}
		
				for (my $n = 0; $n < $K0 / 2; $n++) {
					my $q1 = 2. * PI / $L[1] * $n;
					for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
						next if (($n == 0) && ($m == 0));
	
						my $q2 = 2. * PI / $L[2] * $m;
						my $q = sqrt($q1 * $q1 + $q2 * $q2);
						my $SS = $S4[3 * ($n * ($K1 / 2 + 1) + $m) + 0];
						my $SD = $S4[3 * ($n * ($K1 / 2 + 1) + $m) + 1];
						my $DD = $S4[3 * ($n * ($K1 / 2 + 1) + $m) + 2];
	
						print "$time $q1 $q2 $q $SS $SD $DD\n";
					}
					print "\n";
				}
			}
		}

		if ($vacf_enabled) {
			my $line = <ISF>;
			next unless ($mode == 7);

			chomp $line;
			my @tmp = split / /, $line;
			if ($#tmp < 2) {
				printf("Strange VA: ##$line##\n");
			}
			for (my $k = 0; $k < @tmp; $k++) {
				$va[$i * $max_elements + $j][$k] = $tmp[$k];
			       	$va[$i * $max_elements + $j][$k] /= $c[$i * $max_elements + $j];
			}
		}

		if ($facf_enabled) {
			for (my $k = 0; $k < $facf_rbins; $k++) {
				my $line = <ISF>;
				next unless ($mode == 8);
				chomp $line;
				my @tmp = split / /, $line;

				if ($tmp[2] > 0) {
					$facf_g[$i * $max_elements + $j][$k] = $tmp[0] / $tmp[2];
					$facf_h[$i * $max_elements + $j][$k] = $tmp[1] / $tmp[2];
				} else  {
					$facf_g[$i * $max_elements + $j][$k] = 0;
					$facf_h[$i * $max_elements + $j][$k] = 0;
				}
			}
			for (my $k = 0; $k < $nN; $k++) {
				my $line = <ISF>;
			}
		}

		if ($uacf_enabled) {
			for (my $k = 0; $k < $K0 * ($K1 / 2 + 1); $k++) {
				$line = <ISF>;
				if ($c[$i * $max_elements + $j] > 0 && $mode == 9) {
					$line =~ s/\n//;
					my @tmp = split(/ /, $line);
					$Sr[$k] = $tmp[2];
					$Si[$k] = $tmp[3];

					$Sr[$k] /= $c[$i * $max_elements + $j];
					$Si[$k] /= $c[$i * $max_elements + $j];
				}
			}
	
			if ($c[$i * $max_elements + $j] > 0) {
				if ($mode == 9) {
					print "# t q1 q2 |q| Re Im G\n";
		
					for (my $n = $K0 / 2; $n < $K0; $n++) {
						my $q1 = 2. * PI / $L[1] * ($n - $K0);
						for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
							my $q2 = 2. * PI / $L[2] * $m;
							my $q = sqrt($q1 * $q1 + $q2 * $q2);
							my $re = $Sr[$n * ($K1 / 2 + 1) + $m];
							my $im = $Si[$n * ($K1 / 2 + 1) + $m];
							my $g = $S0[$n * ($K1 / 2 + 1) + $m];

							my $pre = 1. / ($K0 * $K1);
							$re *= $pre * $pre;
							$im *= $pre * $pre;
							$g *= $pre * $pre;

							my $ne = $n - $K0;	
							print "$time $q1 $q2 $q $re $im $g $ne $m\n";
						}
						print "\n";
					}
		
					for (my $n = 0; $n < $K0 / 2; $n++) {
						my $q1 = 2. * PI / $L[1] * $n;
						for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
							my $q2 = 2. * PI / $L[2] * $m;
							my $q = sqrt($q1 * $q1 + $q2 * $q2);
							my $re = $Sr[$n * ($K1 / 2 + 1) + $m];
							my $im = $Si[$n * ($K1 / 2 + 1) + $m];
							my $g = $S0[$n * ($K1 / 2 + 1) + $m];
							
							my $pre = 1. / ($K0 * $K1);
							$re *= $pre * $pre;
							$im *= $pre * $pre;
							$g *= $pre * $pre;
	
							next if (($n == 0) && ($m == 0));
							print "$time $q1 $q2 $q $re $im $g $n $m\n";
						}
						print "\n";
					}
				}
			}
		}

		if ($jacf_enabled && $mode == 11 && $c[$i * $max_elements + $j] > 0) {
			# declare many variables...
 			my @Clsr = (); # l = longitudinal, t = transversal
			my @Cldr = (); # r = real, i = imaginary
			my @Ctsr = (); # s = sum, d = difference
			my @Ctdr = ();
 			my @Clsi = ();
			my @Cldi = ();
			my @Ctsi = ();
			my @Ctdi = ();
			my @q;

			my $pre = $c[$i * $max_elements + $j];

			# read in the stuff
			for (my $n = 0; $n < $K0; $n++) {
				for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
					my $k = $n * ($K1 / 2 + 1) + $m;

					$Clsr[$k] = 0; 
				       	$Clsi[$k] = 0;
					$Cldr[$k] = 0;
				       	$Cldi[$k] = 0; 
					$Ctsr[$k] = 0;
				       	$Ctsi[$k] = 0;
					$Ctdr[$k] = 0;
				       	$Ctdi[$k] = 0;

					$line = <ISF>;
					chomp $line;
					my @tmp = split(/ /, $line);
					splice (@tmp, 0, 8); # drop the FFT values
					if (@tmp != 16) {
						die "JACF: The line '$line' is too short!\n";
					}
					
					next if ($n == 0 && $m == 0);
					$q[0] = ($n < $K0 / 2) ? $n / $L[1] : ($n - $K0) / $L[1];
					$q[1] = $m / $L[2];
					my $inorm = 1. / sqrt($q[0] * $q[0] + $q[1] * $q[1]);
					$q[0] *= $inorm;
				       	$q[1] *= $inorm;
					
					for (my $ti = 0; $ti < 2; $ti++) {
						for (my $tj = 0; $tj < 2; $tj++) {
							my $sr = shift @tmp;
							my $si = shift @tmp;
							my $dr = shift @tmp;
							my $di = shift @tmp;
							
							$sr /= $pre; 
							$si /= $pre; 
							$dr /= $pre; 
							$di /= $pre; 

							$Clsr[$k] += $q[$ti] * $q[$tj] * $sr;
							$Clsi[$k] += $q[$ti] * $q[$tj] * $si;
							$Cldr[$k] += $q[$ti] * $q[$tj] * $dr;
							$Cldi[$k] += $q[$ti] * $q[$tj] * $di;

							$Ctsr[$k] -= $q[$ti] * $q[$tj] * $sr;
							$Ctsi[$k] -= $q[$ti] * $q[$tj] * $si;
							$Ctdr[$k] -= $q[$ti] * $q[$tj] * $dr;
							$Ctdi[$k] -= $q[$ti] * $q[$tj] * $di;

							if ($ti == $tj) {
								$Ctsr[$k] += $sr;
							       	$Ctsi[$k] += $si;
								$Ctdr[$k] += $dr;
							       	$Ctdi[$k] += $di;
							}
						}
					}
				}
			}

			print "# t q1 q2 |q| Cls(t) Cls(0) Cld(t) Cld(0) Cts(t) Cts(0) Ctd(t) Ctd(0) n m\n";
		
			for (my $n = $K0 / 2; $n < $K0; $n++) {
				my $q1 = 2. * PI / $L[1] * ($n - $K0);
				for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
					my $q2 = 2. * PI / $L[2] * $m;
					my $q = sqrt($q1 * $q1 + $q2 * $q2);
					my $k =  $n * ($K1 / 2 + 1) + $m;
					
					$Clsr[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Clsi[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Cldr[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Cldi[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Ctsr[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Ctsi[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Ctdr[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Ctdi[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;

					# the imaginary parts are omitted in the output
					print "$time $q1 $q2 $q ";
					print join (' ',
					       	$Clsr[$k], $Cls0[$k],
					       	$Cldr[$k], $Cld0[$k],
					       	$Ctsr[$k], $Cts0[$k],
					       	$Ctdr[$k], $Ctd0[$k],
						$n - $K0, $m), "\n";
				}
				print "\n";
			}
			for (my $n = 0; $n < $K0 / 2; $n++) {
				my $q1 = 2. * PI / $L[1] * $n;
				for (my $m = 0; $m < $K1 / 2 + 1; $m++) {
					my $q2 = 2. * PI / $L[2] * $m;
					my $q = sqrt($q1 * $q1 + $q2 * $q2);
					my $k =  $n * ($K1 / 2 + 1) + $m;
					next if ($n == 0 && $m == 0);
					
					$Clsr[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Clsi[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Cldr[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Cldi[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Ctsr[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Ctsi[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Ctdr[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;
					$Ctdi[$k] *= ($q1 * $q1 + $q2 * $q2) / $nN;

					# the imaginary parts are omitted in the output
					print "$time $q1 $q2 $q ";
					print join (' ',
						$Clsr[$k], $Cls0[$k],
						$Cldr[$k], $Cld0[$k],
						$Ctsr[$k], $Cts0[$k],
						$Ctdr[$k], $Ctd0[$k],
						$n, $m), "\n";
				}
				print "\n";
			}
		} elsif ($jacf_enabled) {
			# we don't want jacf or the count is zero anyway,
			# thus we just skip this part
			for (my $n = 0; $n < $K0 * ($K1 / 2 + 1); $n++) {
				$line = <ISF>;
			}
		}
	}
}

close(ISF);

if ($mode == 1) { # static collective structure factor
	print "# q1 q2 |q| S\n";

	for (my $i = $K0 / 2; $i < $K0; $i++) {
		my $q1 = 2. * PI / $L[1] * (($i < $K0 / 2) ? $i : $i - $K0);
		for (my $j = 0; $j < $K1 / 2 + 1; $j++) {
			my $q2 = 2. * PI / $L[2] * $j;
			my $q = sqrt($q1 * $q1 + $q2 * $q2);
			my $S = $G[$i * ($K1 / 2 + 1) + $j];
			print "$q1 $q2 $q $S\n";
		}
		print "\n";
	}

	for (my $i = 0; $i < $K0 / 2; $i++) {
		my $q1 = 2. * PI / $L[1] * (($i < $K0 / 2) ? $i : $i - $K0);
		for (my $j = 0; $j < $K1 / 2 + 1; $j++) {
			my $q2 = 2. * PI / $L[2] * $j;
			my $q = sqrt($q1 * $q1 + $q2 * $q2);
			my $S = $G[$i * ($K1 / 2 + 1) + $j];
			print "$q1 $q2 $q $S\n";
		}
		print "\n";
	}
}

if ($mode == 4) { # mean-square displacement
	print "# t D alpha\n";

	for (my $i = 0; $i < $blocks; $i++) {
		my $len = min($blocklen[$i], $max_elements);
		for (my $j = 0; $j < $len; $j++) {
			my $time = ($j + 1) * $dt * (($max_elements + 1) ** $i);
			my $count = $c[$i * $max_elements + $j];

			if ($count > 0) {
				my $ir2 = $r2[$i * $max_elements + $j];
				$ir2 /= $nN * $count;
				my $ir4 = $r4[$i * $max_elements + $j];
				$ir4 /= $nN * $count;

				my $a = .5 * $ir4 / ($ir2 * $ir2) - 1.;
				print "$time $ir2 $a\n";
			}
		}
	}
}

if ($mode == 5) { # Q_4 and chi-4
	print "# t <Q_S> <Q_D> <Q_S*Q_S> <Q_S*Q_D> <Q_D*Q_D>\n";
	
	for (my $i = 0; $i < $blocks; $i++) {
		my $len = min($blocklen[$i], $max_elements);
		for (my $j = 0; $j < $len; $j++) {
			my $time = ($j + 1) * $dt * (($max_elements + 1) ** $i);
			my $offset = $i * $max_elements + $j;
			my $count = $c[$offset];

			if ($count > 0) {
				print "$time ";
				for (my $k = 0; $k < @{$fp_Q[$offset]}; $k++) {
					my $q = $fp_Q[$offset][$k];
					$q /= $count;
					print "$q ";
				}
				print "\n";
			}
		}
	}
}

if ($mode == 7) { # velocity autocorrelation function
	print "# t <v_x(t)v_x(0)> <v_y(t)v_y(0)> <v_z(t)v_z(0)>\n";

	for (my $i = 0; $i < $blocks; $i++) {
		my $len = min($blocklen[$i], $max_elements);
		for (my $j = 0; $j < $len; $j++) {
			my $time = ($j + 1) * $dt * (($max_elements + 1) ** $i);
			my $offset = $i * $max_elements + $j;
			my $count = $c[$offset];
			
			if ($count > 0) {
				print "$time ";
				for (my $k = 0; $k < @{$va[$offset]}; $k++) {
					my $v = $va[$offset][$k];
					$v /= $nN;
					print "$v ";
				}
				print "\n";
			}
		}
	}
}

if ($mode == 8) { # force autocorrelation function
	print "# t r g(t,r) h(t,r)\n";

	for (my $i = 0; $i < $blocks; $i++) {
		my $len = min($blocklen[$i], $max_elements);
		for (my $j = 0; $j < $len; $j++) {
			my $time = ($j + 1) * $dt * (($max_elements + 1) ** $i);
			my $offset = $i * $max_elements + $j;
			my $count = $c[$offset];
			
			if ($count > 0) {
				for (my $k = 0; $k < $facf_rbins; $k++) {
					my $gtr = $facf_g[$offset][$k];
					my $htr = $facf_h[$offset][$k];
					my $r = $k * $facf_dr;
					print "$time $r $gtr $htr\n";
				}
				print "\n";
			}
		}
	}
}

if ($mode == 10) { # static undulation power spectrum
	print "# q1 q2 |q| S\n";

	for (my $i = $K0 / 2; $i < $K0; $i++) {
		my $q1 = 2. * PI / $L[1] * (($i < $K0 / 2) ? $i : $i - $K0);
		for (my $j = 0; $j < $K1 / 2 + 1; $j++) {
			my $q2 = 2. * PI / $L[2] * $j;
			my $q = sqrt($q1 * $q1 + $q2 * $q2);
			my $S = $S0[$i * ($K1 / 2 + 1) + $j];
			$S *= 1. / ($K0 * $K1);
			$S *= 1. / ($K0 * $K1);
			print "$q1 $q2 $q $S\n";
		}
		print "\n";
	}

	for (my $i = 0; $i < $K0 / 2; $i++) {
		my $q1 = 2. * PI / $L[1] * (($i < $K0 / 2) ? $i : $i - $K0);
		for (my $j = 0; $j < $K1 / 2 + 1; $j++) {
			my $q2 = 2. * PI / $L[2] * $j;
			my $q = sqrt($q1 * $q1 + $q2 * $q2);
			my $S = $S0[$i * ($K1 / 2 + 1) + $j];
			$S *= 1. / ($K0 * $K1);
			$S *= 1. / ($K0 * $K1);
			print "$q1 $q2 $q $S\n";
		}
		print "\n";
	}
}


if ($mode == 12) { # static current power spectrum
	print "# q1 q2 |q| Cls Cld Cts Ctd\n";

	for (my $i = $K0 / 2; $i < $K0; $i++) {
		my $q1 = 2. * PI / $L[1] * (($i < $K0 / 2) ? $i : $i - $K0);
		for (my $j = 0; $j < $K1 / 2 + 1; $j++) {
			my $q2 = 2. * PI / $L[2] * $j;
			my $q = sqrt($q1 * $q1 + $q2 * $q2);
			my $k = $i * ($K1 / 2 + 1) + $j;

			my $ls = $Cls0[$k];
			my $ld = $Cld0[$k];
			my $ts = $Cts0[$k];
			my $td = $Ctd0[$k];
			print "$q1 $q2 $q $ls $ld $ts $td\n";
		}
		print "\n";
	}

	for (my $i = 0; $i < $K0 / 2; $i++) {
		my $q1 = 2. * PI / $L[1] * (($i < $K0 / 2) ? $i : $i - $K0);
		for (my $j = 0; $j < $K1 / 2 + 1; $j++) {
			my $q2 = 2. * PI / $L[2] * $j;
			my $q = sqrt($q1 * $q1 + $q2 * $q2);
			my $k = $i * ($K1 / 2 + 1) + $j;
			next if ($i == 0 && $j == 0);

			my $ls = $Cls0[$k];
			my $ld = $Cld0[$k];
			my $ts = $Cts0[$k];
			my $td = $Ctd0[$k];
			print "$q1 $q2 $q $ls $ld $ts $td\n";
		}
		print "\n";
	}
}

