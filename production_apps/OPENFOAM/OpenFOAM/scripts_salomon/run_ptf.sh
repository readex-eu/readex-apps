#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N FOAM_PTF
#PBS -l select=2:ncpus=24:mpiprocs=24:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=03:00:00
#PBS -m be

APP="simpleFoam -parallel"
THRDS=1
MPI_PROCS=24
PHASE_REG_NAME="iteration"

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
export FM_DIR=$FM_DIR/..

cd $FM_DIR
source readex_env/set_env_ptf_rapl.source
source scripts_$READEX_MACHINE/environment.sh

cd $FM_DIR/OpenFOAM-v1612+/
source $FM_DIR/OpenFOAM-v1612+/etc/bashrc
cd ../../motorBike24/

export SCOREP_TOTAL_MEMORY=3G

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin
export SCOREP_RRL_VERBOSE="WARN"

export SCOREP_METRIC_PLUGINS=x86_energy_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN=*/E
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_MPI_ENABLE_GROUPS=ENV

psc_frontend --apprun="$APP" --mpinumprocs=$MPI_PROCS --ompnumthreads=$THRDS --phase=$PHASE_REG_NAME --tune=readex_intraphase --config-file=readex_config_extended.xml --force-localhost --info=2 --selective-info=AutotuneAll,AutotunePlugins

