#!/bin/bash

name=$1

files="viscosity_0?_$name/viscosity.*.dat"
prog=~/fdmdpd2/src/viscosity/int.pl
end=$2
dt=0.005

for filename in $files; do
	echo -n "$filename " >> er_$name
	~/fdmdpd2/src/viscosity/vis 100000 500 8 $filename | $prog $dt 0 $end 1 >> er_$name
done


