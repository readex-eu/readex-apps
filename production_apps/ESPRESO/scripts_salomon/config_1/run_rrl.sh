#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N espreso_rrl
#PBS -l select=1:ncpus=24:mpiprocs=2:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=24:00:00  
#PBS -m be

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
cd $FM_DIR/../../

source env/readex_env/set_env_rrl.source
source env/modules.taurus

source env/paths.default
source env/threading.default 12

export SCOREP_SUBSTRATE_PLUGINS='rrl'
export SCOREP_RRL_VERBOSE="WARN"
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin,OpenMPTP
export SCOREP_RRL_TMM_PATH=tuning_model.json
export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=false
export SCOREP_MPI_ENABLE_GROUPS=ENV

echo "running for RRL"
mpirun -n 2 ./espreso_rdd -c espreso.ecf

echo "running for RRL done"
