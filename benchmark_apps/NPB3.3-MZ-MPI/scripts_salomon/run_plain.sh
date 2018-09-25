#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N BTMZ_plain
#PBS -o BTMZ_plain.out
#PBS -e BTMZ_plain.err
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
source readex_env/set_env_plain.source
module list

# run application
mpirun $BIND_TO_SOCKETS ./bin/bt-mz.C.2_plain


