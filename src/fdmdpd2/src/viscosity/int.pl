#!/usr/bin/perl -w
use strict;

sub syntax
{
	print STDERR "Syntax: $0 [h] [a] [b] [column]\n";
	print STDERR "h: binwidth, a: lower boundary, b: upper boundary\n";
	print STDERR "column: column with function values\n";
	print STDERR "a and b start at 0.\n\n";
	exit 1;
}

my $h = shift; # timestep
my $a = shift; # first bin
my $b = shift; # last bin
my $col = shift; # column

syntax() if (!defined $h || !defined $a || !defined $b || !defined $col);

my $int = 0.;
my @f;

while(<>) {
	next if /^#/;
	my @tmp = split / /;
	push @f, $tmp[$col];
}

#print "\@f hat ". @f ." Elemente\n";
splice(@f,0,$a);
splice(@f,$b-$a);
#print "\@f hat jetzt noch ". @f ." Elemente\n";

$int =  shift(@f) * 3. / 8.;
$int += shift(@f) * 7. / 6.;
$int += shift(@f) * 23. / 24.;

$int += pop(@f) * 3. / 8.;
$int += pop(@f) * 7. / 6.;
$int += pop(@f) * 23. / 24.;

foreach my $k (@f) {
	$int += $k;
}

print $int * $h, "\n";
