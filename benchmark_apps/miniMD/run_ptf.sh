#!/bin/sh

#SBATCH --time=5:00:00   # walltime
#SBATCH --nodes=2  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -J "e_ptf"   # job name
#SBATCH -A p_readex

echo "run PTF begin."

NP=1 # check against --ntasks and tasks-per-node

module purge
source ./readex_env/set_env_ptf_hdeem.source

INPUT_FILE=in3.data #in.lj.miniMD
PHASE=INTEGRATE_RUN_LOOP

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

#export SCOREP_METRIC_PLUGINS=x86_energy_sync_plugin
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN="*/E"
#export SCOREP_METRIC_PLUGINS_SEP=";"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_CONNECTION="INBAND"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE="WARN"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_MPI_ENABLE_GROUPS=ENV

psc_frontend --apprun="./miniMD_openmpi_ptf -i $INPUT_FILE" --mpinumprocs=$NP --ompnumthreads=24 --phase=$PHASE --tune=readex_intraphase --config-file=ptf_readex_config.xml --force-localhost --info=2 --selective-info=AutotuneAll,AutotunePlugins

echo "run PTF done."
