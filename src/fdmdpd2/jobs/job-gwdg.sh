#!/bin/bash
# Jobscript for the Nehalem Cluster at GWDG
#
#BSUB -J jobname
#BSUB -q gwdg-x64par2
#BSUB -a intelmpi
#BSUB -R "model==IE5540"
#BSUB -wa 'URG'
#BSUB -wt '2'
#BSUB -u mhoembe@gwdg.de
#BSUB -N
#BSUB -o stdout.%J.%I
#BSUB -e stderr.%J.%I
#BSUB -W 01:00
#BSUB -n 8

export LD_LIBRARY_PATH=/opt/intel/Compiler/11.1/059/lib/intel64 
mpirun.lsf ./fdmdpd2 4 2 config.txt resume.dat
