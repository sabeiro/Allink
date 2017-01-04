#!/bin/bash

files="corr.out.*"

for file in $files; do
	cat $file | ~/fdmdpd2/src/viscosity/int.pl 0.005 0 4000 1
done

