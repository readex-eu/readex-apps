#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N espreso_ptf
#PBS -l select=2:ncpus=24:mpiprocs=2:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=24:00:00  
#PBS -m be

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
cd $FM_DIR/../../

source env/readex_env/set_env_ptf_rapl.source
source env/environment.sh

source env/paths.default
source env/threading.default 12

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin,OpenMPTP
export SCOREP_RRL_VERBOSE="WARN"

#export SCOREP_METRIC_PLUGINS=hdeem_sync_plugin
#export SCOREP_METRIC_PLUGINS_SEP=";"
#export SCOREP_METRIC_HDEEM_SYNC_PLUGIN="*/E"
#export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_CONNECTION="INBAND"
#export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_VERBOSE="WARN"
#export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_METRIC_PLUGINS=x86_energy_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN=*/E
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_MPI_ENABLE_GROUPS=ENV

PHASE=Main

psc_frontend --apprun="./espreso_ptf -c espreso.ecf" --mpinumprocs=2 --ompnumthreads=12 --phase=$PHASE --tune=readex_intraphase --config-file=readex_config_extended.xml --force-localhost --info=1 --selective-info=AutotuneAll,AutotunePlugins,FrontendStateMachines,AutotuneAgentStrategy,ApplTuningParameter | tee -a LOGptf


