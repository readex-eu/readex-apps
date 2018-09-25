#!/bin/sh

#SBATCH --time=24:00:00   # walltime
#SBATCH --nodes=5  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=READEX
#SBATCH -J "amg2013_ptf"   # job name
#SBATCH --output=amg2013_ptf.out
#SBATCH --error=amg2013_ptf.out

cd ..

module purge
source ./readex_env/set_env_ptf_hdeem.source
#source ./readex_env/set_env_ptf_rapl.source

export PSC_CPU_BIND="--cpu_bind=verbose,sockets"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin,OpenMPTP
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

export OMP_NUM_THREADS=12

PHASE=main_phase

psc_frontend --apprun="./test/amg2013_ptf -P 2 2 2 -r 40 40 40" --mpinumprocs=8 --ompnumthreads=12 --phase=$PHASE --tune=readex_intraphase --config-file=readex_config_ptf.xml --force-localhost --info=2 --selective-info=AutotuneAll,AutotunePlugins

