#!/bin/bash
#PBS -N fep
#PBS -M mhoembe@gwdg.de
#PBS -A nic00015
#PBS -m e
#PBS -l naccesspolicy=singleuser
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=4
#PBS -l feature=xe

BINDIR=$HOME/fdmdpd2/src
FDMDPD2=$BINDIR/fdmdpd2
O2N=$BINDIR/scripts/o2n.pl
np=$(cat $PBS_NODEFILE | wc -l)

cd $WORK/fdmdpd2/RAFT-DIFFUSION/rdiff-100-30-14
#list=output000[01234]*.dat
#list="output000[56789]*.dat"
#list="output001[01234]*.dat"
list="output001[56789]*.dat"

for file in $list; do
	$O2N $file > tmpnew.dat
	mpiexec -v -np $np $FDMDPD2 2 2 ../config-fep.dat tmpnew.dat $HOME/fdmdpd2/nvt ../virials.txt 
	echo $file
done
rm -f tmpnew.dat
exit 0
