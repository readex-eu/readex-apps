#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_ptf
#PBS -l select=2:ncpus=24:mpiprocs=24:accelerator=false
#PBS -l walltime=05:00:00

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
    LOC=$PBS_O_WORKDIR
else
    LOC=$(pwd)
fi

cd $LOC

module load mkl
module load CMake/3.9.1

. readex_env/set_env_ptf_rapl.source
. ../paths.source

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin
export SCOREP_RRL_VERBOSE="WARN"
export SCOREP_MPI_ENABLE_GROUPS=ENV

export SCOREP_METRIC_PLUGINS=x86_energy_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN=*/E
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000


################################################################################

cd $LOC/../ELMER-APP/

NPROCS=4

psc_frontend --apprun="ElmerSolver_mpi" --mpinumprocs=${NPROCS} --ompnumthreads=1 --phase="ElmerSolver_core" --info=2 --config-file=readex_config_extended.xml --tune=readex_config --readex-app-config="readex_acp_config.cfg"  --force-localhost --selective-info=AutotunePlugins,FrontendStateMachines 2>&1 | tee -a log_JOB.txt




