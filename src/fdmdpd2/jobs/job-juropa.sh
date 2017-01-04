#!/bin/bash
#PBS -N rho40_kN100_chiN30-nvt
#PBS -M mhoembe@gwdg.de
#PBS -m e
#PBS -l naccesspolicy=singlejob
#PBS -l walltime=6:00:00
#PBS -l nodes=8:ppn=8
#PBS -l signal=15@120

cd $PBS_O_WORKDIR
np=$(cat $PBS_NODEFILE | wc -l)
mpiexec -v -np $np ./fdmdpd2 8 8 config.txt resume.dat
