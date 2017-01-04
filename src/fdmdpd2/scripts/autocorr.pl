#!/usr/bin/perl -w

# read data into @x, skipping infinities
my @x;
while (<>) {
    if ($_ < -1e99 || $_ > +1e99) {
	if (@x) {
	    push @x, $x[$#x];
	}
    } else {
	push @x, $_;
    }
}

# compute mean, variance
my $n = @x;
my ($sum, $sqsum) = (0, 0);
grep ($sum += $_, @x);
grep ($sqsum += $_**2, @x);
my $mean = $sum / $n;
my $variance = $sqsum / $n - $mean**2;

# print autocorrelation
for my $lag (0 .. $n-1) {
    my ($prod, $count) = (0, 0);
    for my $i (0 .. $n-1-$lag) {
	$prod += ($x[$i] - $mean) * ($x[$i+$lag] - $mean);
	++$count;
    }
    print $lag, " ", $prod/$count/$variance, "\n";
}

