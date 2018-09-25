#!/bin/bash

#SBATCH --time=02:30:00
#SBATCH --nodes=2
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH -p haswell
#SBATCH --mem-per-cpu=600M   # memory per CPU core
#SBATCH -J "READEX_kripke"
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH --comment="no_monitoring"
#SBATCH --mail-type=ALL

cd ../build
. ../readex_env/set_env_ptf_hdeem.source
. ../environment.sh

#cp ../RESULTS/readex_config_extended.xml .

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin
export SCOREP_RRL_VERBOSE="WARN"

export SCOREP_METRIC_PLUGINS=hdeem_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN="*/E"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000
clearHdeem

#export SCOREP_METRIC_PLUGINS=x86_energy_sync_plugin
#export SCOREP_METRIC_PLUGINS_SEP=";"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN=*/E
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_CONNECTION="INBAND"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE="WARN"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_MPI_ENABLE_GROUPS=ENV

# run the application
psc_frontend --apprun="./kripke $KRIPKE_COMMAND" --mpinumprocs=24 --ompnumthreads=1 --phase="Loop" --info=2 --config-file=readex_config_extended.xml --tune=readex_intraphase  --force-localhost --selective-info=AutotunePlugins,FrontendStateMachines 2>&1 | tee -a log_JOB.txt

echo "running psc_frontend done"

