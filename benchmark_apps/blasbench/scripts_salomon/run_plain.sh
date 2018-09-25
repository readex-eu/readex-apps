#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N blasbench_plain
#PBS -o blasbench_plain.out
#PBS -e blasbench_plain.err
#PBS -l select=1:ncpus=24:mpiprocs=2:ompthreads=12:accelerator=true
#PBS -l walltime=00:30:00
###############################################################################

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
  LOC=$PBS_O_WORKDIR
else
  LOC=$(pwd)
fi

cd $LOC/..

module purge
source ./readex_env/set_env_plain.source
source $LOC/set_env_blasbench.source

mpirun $BIND_TO_SOCKETS ./blasbench_plain --dgemv 300,8192 --dgemm 150,2048 --tan 60000 --phase 10 

