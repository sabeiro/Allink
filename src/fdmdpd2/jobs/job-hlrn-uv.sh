#!/bin/bash
#PBS -N jobname
#PBS -M mhoembe@gwdg.de
#PBS -m e
#PBS -l walltime=12:00:00
#PBS -l nodes=4:ppn=8
#PBS -l signal=23@300
#PBS -l feature=uv

. ${MODULESHOME}/init/bash
module load mpt

cd $PBS_O_WORKDIR
mpiexec -verbose -np 32 ../fdmdpd2 8 4 config.txt resume.dat

