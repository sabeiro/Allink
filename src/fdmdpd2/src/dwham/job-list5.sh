#!/bin/bash
#PBS -N dwham-juropa
#PBS -M mhoembe@gwdg.de
#PBS -m e
#PBS -l naccesspolicy=singlejob
#PBS -l walltime=6:00:00
#PBS -l nodes=128:ppn=8

export MV2_CPU_MAPPING="0:1:4:5:2:3:6:7"
export MV2_ENABLE_AFFINITY="1"
export MV2_NUM_PORTS="2"
export LD_LIBRARY_PATH="$HOME/local/lib:$LD_LIBRARY_PATH"

cd $PBS_O_WORKDIR
np=$(cat $PBS_NODEFILE | wc -l)
mkdir $PBS_JOBNAME-$PBS_JOBID	
cd $PBS_JOBNAME-$PBS_JOBID

cp ../f-list5 ./.f.dat

mpiexec -v -np $np ../dwham ../config.dat ../list5

if [ "$?" -ne 0 ]; then
        echo "Job failed."
        exit 1
fi

echo "Job finished."
exit 0
