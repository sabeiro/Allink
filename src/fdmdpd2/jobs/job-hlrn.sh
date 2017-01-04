#!/bin/bash
# This is a sample job script for FDMDPD2 at HLRN
# 
#PBS -N jobname
#PBS -M user@theorie.physik.uni-goettingen.de
#PBS -A nic00016
#PBS -m e
#PBS -l walltime=24:00:00
#PBS -l nodes=3:ppn=8
#PBS -l feature=xe
#PBS -l signal=15@120

JOB="job-hlrn.sh"

resubmit()
{
	if [ -f $JOB ]; then
		msub $JOB
	fi
}

on_term()
{
	echo "Script caught signal. Resubmitting..."
	resubmit
	exit 0
}

trap 'on_term' USR1
trap 'on_term' TERM

cd $PBS_O_WORKDIR
np=$(cat $PBS_NODEFILE | wc -l)
mpiexec -v -np $np ../fdmdpd2 3 8 config.dat resume.dat

if [ "$?" -ne 0 ]; then
        echo "Job failed"
        exit 1
fi

echo "Script finished normally. Resubmitting..."
resubmit
exit 0
