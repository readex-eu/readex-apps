#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N espreso_plain
#PBS -l select=1:ncpus=24:mpiprocs=2:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=24:00:00  
#PBS -m be

################################################################################
if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
cd $FM_DIR/../../

source env/readex_env/set_env_plain.source
source env/modules.taurus

source env/paths.default
source env/threading.default 12

# first run without measurement
mpirun -n 2 --exclusive -p interactive --mem-per-cpu 2500M ./espreso_plain -c espreso.ecf

