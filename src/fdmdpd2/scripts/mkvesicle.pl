#!/usr/bin/perl -w
#
# mkvesicle.pl - Create a vesicle initial configuration
#
# The algorithm to distribute the lipids over the sphere is called
# the "Greek spiral". See http://www.xsi-blog.com/archives/115
#
use strict;

my $NA = 8;
my $NB = 2;
my $Re = 2.0;
my $kb = 3.;
my $ks = 20.537;
my $rho = 5.0;
my $kN= 70.0;
my $chiN = 30.0;
my $VIRIAL = "~/fdmdpd2/scripts/calc-virials.pl";

my $ri = shift or syntax();
my $ni = shift or syntax();
my $ra = shift or syntax();
my $na = shift or syntax();

my @lipids = ();

sub syntax
{
	print STDERR "$0 [R_i] [n_i] [R_a] [n_a]\n";
	print STDERR "R_i / R_a: Radius of the inside / outside\n";
	print STDERR "n_i / n_a: # lipids on the inner / outer shell\n\n";
	exit 1;
}

# print header
my $l = 5. * $ra;
my $n = $ni + $na;
my $N = $NA + $NB;
print "# L=$l $l $l t=0 blocks=2\n";
print `$VIRIAL $rho $kN $chiN`;
print "# a2=0.9 a3=1 Re=$Re N=$N ks=$ks kb=$kb l0=0\n";
print "# n=$ni N=$N name=INSIDE\n";

# put lipids on the inside
my $inc = 3.141592654 * (3. - sqrt(5.));
my $off = 2. / $ni;
for (my $i = 0; $i < $ni; $i++) {
	my $y = $i * $off - 1. + ($off / 2.);
	my $r = sqrt(1. - $y * $y);
	my $phi = $i * $inc;

	my $x = cos($phi) * $r;
	my $z = sin($phi) * $r;

	put_lipids(-1., $x*$ri, $y*$ri, $z*$ri);
}
print_lipids();
@lipids = ();

# put lipids on the outside
print "# n=$na N=$N name=OUTSIDE\n";
$off = 2. / $na;
for (my $i = 0; $i < $na; $i++) {
	my $y = $i * $off - 1. + ($off / 2.);
	my $r = sqrt(1. - $y * $y);
	my $phi = $i * $inc;

	my $x = cos($phi) * $r;
	my $z = sin($phi) * $r;

	put_lipids(1., $x*$ra, $y*$ra, $z*$ra);
}
print_lipids();

sub gauss
{
	my $u1 = rand();
	my $u2 = rand();
	sqrt(-2. * log($u1)) * cos(2. * 3.141592654 * $u2);
}

sub print_lipids
{
	for (my $i = 0; $i < @lipids; $i++) {
		my $vx = gauss();
		my $vy = gauss();
		my $vz = gauss();
		print join(' ', $lipids[$i]{'x'}, $lipids[$i]{'y'}, $lipids[$i]{'z'}, $vx, $vy, $vz, $lipids[$i]{'t'}), "\n";
	}
}

sub put_lipids
{
	my $pref = shift;
	my $x = shift;
	my $y = shift;
	my $z = shift;
	my @lipid = ();

	# last A bead goes onto the sphere
	$lipid[$NA - 1]->{'x'} = $x;
	$lipid[$NA - 1]->{'y'} = $y;
	$lipid[$NA - 1]->{'z'} = $z;
	$lipid[$NA - 1]->{'t'} = 0;

	my $phi = atan2($y, $x);
	my $r = sqrt($x*$x + $y*$y + $z*$z);
	my $ct = $z / $r;

	# now put the B beads
	for (my $i = $NA; $i < $NA + $NB; $i++) {
		$phi += (rand() - .5) * 2. * 3.141592654 * 0.01;
		$r += 0.2 * $pref;
		$lipid[$i]->{'x'} = $r * cos($phi) * sqrt(1. - $ct*$ct);
		$lipid[$i]->{'y'} = $r * sin($phi) * sqrt(1. - $ct*$ct);
		$lipid[$i]->{'z'} = $r * $ct;
		$lipid[$i]->{'t'} = 1;
	}
	
	$phi = atan2($y, $x);
	$r = sqrt($x*$x + $y*$y + $z*$z);

	# other A beads
	for (my $i = 0; $i < $NA - 1; $i++) {
		$phi += (rand() - .5) * 2. * 3.141592654 * 0.01;
		$r -= 0.2 * $pref;
		$lipid[$i]->{'x'} = $r * cos($phi) * sqrt(1. - $ct*$ct);
		$lipid[$i]->{'y'} = $r * sin($phi) * sqrt(1. - $ct*$ct);
		$lipid[$i]->{'z'} = $r * $ct;
		$lipid[$i]->{'t'} = 0;
	}

	for (my $i = 0; $i < $NA + $NB; $i++) {
		push @lipids, $lipid[$i];
	}
}

