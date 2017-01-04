#!/bin/bash
# This is a sample job script for FDMDPD2 at the Rocks-Cluster
# 
#PBS -N fdmdpd2
#PBS -M mhoembe@gwdg.de
#PBS -m e
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=8
#PBS -l vmem=8gb
#PBS -l mem=8gb

cd $PBS_O_WORKDIR
np=$(cat $PBS_NODEFILE | wc -l)

if [ -e .lock ]; then
	echo ".lock exists -- previous job didn't clean up properly."
	exit 1
else
	date > .lock
	env >> .lock
	cat $PBS_NODEFILE >> .lock

	# the job terminates cleanly after 86000 seconds (approx. 23:54 hours)
	echo "sleep 86000" > killme
	echo "qsig -s 23 $PBS_JOBID" >> killme
	chmod 755 killme
	./killme &
	
	mpiexec -v -np $np ../fdmdpd2 4 2 config.txt resume.dat
	if [ "$?" != "0" ]; then
		echo "FDMDPD2 failed."
		exit 1
	else
		echo "FDMDPD2 succeeded."
		rm -f .lock
		if [ -e job.sh ]; then
			qsub job.sh
		fi
	fi
	exit 0
fi

