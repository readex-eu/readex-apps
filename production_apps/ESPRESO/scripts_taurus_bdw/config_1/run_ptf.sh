#!/bin/sh

#SBATCH -t 1-12:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ondrej.vysocky@vsb.cz

#SBATCH --nodes=2
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=14

#SBATCH --partition=broadwell 
#SBATCH --reservation=p_readex_56
#SBATCH -A p_readex	#to account your compute time on the readex project
#SBATCH --exclusive	
#SBATCH --mem-per-cpu=2200M
#SBATCH --comment="no_monitoring"

cd ../../

source env/readex_env/set_env_ptf.source
source env/environment.sh
source env/paths.default
source env/threading.default 14

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

psc_frontend --apprun="./espreso_ptf -c espreso.ecf" --mpinumprocs=2 --ompnumthreads=14 --phase=$PHASE --tune=readex_intraphase --config-file=readex_config_extended.xml --force-localhost --info=1 --selective-info=AutotuneAll,AutotunePlugins,FrontendStateMachines,AutotuneAgentStrategy,ApplTuningParameter


