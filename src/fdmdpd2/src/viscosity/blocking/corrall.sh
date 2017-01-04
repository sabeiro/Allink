#!/bin/bash

files="out.*"

for file in $files; do
#	cat $file | ~/fdmdpd2/src/corr/corr -C -D 0.005 -L 50000 -A 8 -B 8 2>/dev/null | ~/fdmdpd2/src/viscosity/int.pl 0.005 0 10000 1
	cat $file | ~/fdmdpd2/src/corr/corr -C -D 0.005 -L 50000 -A 8 -B 8 2>/dev/null > corr.$file
done

