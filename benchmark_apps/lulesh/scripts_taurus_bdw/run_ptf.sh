#!/bin/sh

#SBATCH --time=0:30:00   # walltime

#SBATCH --nodes=2  # number of processor cores (i.e. tasks)
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=28

#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --mem-per-cpu=2000M   # memory per CPU core
#SBATCH --reservation=p_readex_56
#SBATCH -J "e_ptf_bdw"   # job name
#SBATCH -A p_readex
#SBATCH --output=lulesh_ptf.out
#SBATCH --error=lulesh_ptf.out
#SBATCH --comment="no_monitoring"

cd ..
echo "run PTF begin."

module purge
source ./readex_env/set_env_ptf.source 

PHASE=foo

export LD_LIBRARY_PATH+=":/usr/local/lib/"

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin,OpenMPTP"
export SCOREP_RRL_VERBOSE="WARN"
export SCOREP_TUNING_CPU_FREQ_PLUGIN_VERBOSE=DEBUG
export SCOREP_TUNING_UNCORE_FREQ_PLUGIN_VERBOSE=DEBUG
export SCOREP_TUNING_OPENMPTP_PLUGIN_VERBOSE=DEBUG

export SCOREP_METRIC_PLUGINS=x86_energy_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN="*/E"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

#export SCOREP_METRIC_PLUGINS=x86_energy_sync_plugin
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN="*/E"
#export SCOREP_METRIC_PLUGINS_SEP=";"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_CONNECTION="INBAND"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE="WARN"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000
export OMP_NUM_THREADS=28
#export GOMP_CPU_AFFINITY=0-23
export SCOREP_TOTAL_MEMORY=3G
export SCOREP_MPI_ENABLE_GROUPS=ENV

psc_frontend --apprun="./lulesh2.0_ptf -s 150" --mpinumprocs=1 --ompnumthreads=28 --phase=$PHASE --tune=readex_intraphase --config-file=readex_config_ptf.xml --force-localhost --info=2 --selective-info=AutotuneAll,AutotunePlugins

echo "run PTF done."
