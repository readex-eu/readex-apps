#!/bin/sh

#SBATCH --time=6:00:00   # walltime
#SBATCH --nodes=2  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=READEX
#SBATCH -J "blasbench_ptf"   # job name
#SBATCH --output=blasbench_ptf.out
#SBATCH --error=blasbench_ptf.out

LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_ptf_hdeem.source
#source ./readex_env/set_env_ptf_rapl.source
source $LOC/set_env_blasbench.source

export PSC_CPU_BIND="--cpu_bind=verbose,sockets"

export SCOREP_TOTAL_MEMORY=2G

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

PHASE=Main

psc_frontend --apprun="./blasbench_ptf --dgemv 300,8192 --dgemm 150,2048 --tan 60000 --phase 75" --mpinumprocs=2 --ompnumthreads=12 --phase=$PHASE --tune=readex_intraphase --config-file=readex_config_ptf.xml --force-localhost --info=2 --selective-info=AutotuneAll,AutotunePlugins

