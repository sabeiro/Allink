#!/usr/bin/perl -w
# 
# separate-species.pl - Move different molecular species into their own blocks
#
use strict;

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
my $str;
my $arch;
my $otherheaders = '';
my $beadcnt;
my @lines;
my @archcnt;

while(<>) {
	if (/^# L=([^ ]+) ([^ ]+) ([^ ]+) t=([^ ]+) blocks=(\d+)/) {
		print STDERR "$_";
		$L[0] = $1;
		$L[1] = $2;
		$L[2] = $3;
		$time = $4;
		$blocks = $5;
		$nN = 0;
		$block = -1;
		$c = 0;
		for (my $i = 0; $i < $blocks; $i++) {
			foreach my $a (keys %{$lines[$i]}) {
				$lines[$i]{$a} = '';
				$archcnt[$i]{$a} = 0;
			}
		}
		$otherheaders = '';
	}
	elsif (/^# n=(\d+) N=(\d+) name=(\w+)/) {
		$block++;
		$nN += $1 * $2;
		$n[$block] = $1;
		$N[$block] = $2;
		$resname[$block] = $3;
		$beadcnt = 0;
	}
	elsif (/^#/) { 
		$otherheaders .= $_;
	}
	else {
		my $l = $_;
		$str .= $l;
		chomp $l;
		my @tmp = split / /, $l;
		$arch .= $tmp[6];
		$beadcnt++;
		$c++;

		if ($beadcnt == $N[$block]) {
#			print "$arch\n";
			$lines[$block]{$arch} .= $str;
			$archcnt[$block]{$arch}++;
			$str = '';
			$beadcnt = 0;
			$arch = '';
		}
	}

	if ($c == $nN && $block + 1 == $blocks) {
		my $newbl = 0;
		for (my $i = 0; $i < $blocks; $i++) {
			$newbl += keys %{$lines[$i]};
		}

		print "# L=$L[0] $L[1] $L[2] t=$time blocks=$newbl\n";
		print $otherheaders;
		for (my $i = 0; $i < $blocks; $i++) {
			my $r = 0;
			foreach my $a (sort keys %{$lines[$i]}) {
				print "# n=$archcnt[$i]{$a} N=$N[$i] name=$resname[$i]_$r\n";
				print $lines[$i]{$a};
				$r++;
			}
		}
	}
}

