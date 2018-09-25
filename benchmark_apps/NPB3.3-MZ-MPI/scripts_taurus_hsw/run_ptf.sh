#!/bin/sh

#SBATCH --time=08:00:00   # walltime
#SBATCH --nodes=2 # number of processor cores (i.e. tasks)
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --reservation=READEX
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "NPB_C_dta"   # job name
#SBATCH --output=bt_dta.out
#SBATCH --error=bt_dta.out
###############################################################################

ps aux | grep diamon

cd ..
module purge
source readex_env/set_env_ptf_hdeem.source
module list

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin,OpenMPTP"
export SCOREP_RRL_VERBOSE="WARN"

export SCOREP_METRIC_PLUGINS=hdeem_sync_plugin
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN="*/E"
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true

export SCOREP_TUNING_CPU_FREQ_PLUGIN_VERBOSE=DEBUG
export SCOREP_TUNING_UNCORE_FREQ_PLUGIN_VERBOSE=DEBUG
export OMP_NUM_THREADS=12
PHASE=phase
app=bt-mz.C.24_ptf
psc_frontend --apprun="./bin/bt-mz.C.2_ptf" \
             --mpinumprocs=2 \
             --ompnumthreads=12 \
             --phase=$PHASE \
             --tune=readex_intraphase \
             --config-file=readex_config.xml \
             --force-localhost \
             --info=7 \
             --selective-info=AutotuneAll,AutotunePlugins,FrontendStateMachines,AutotuneAgentStrategy \

