#!/bin/sh

#SBATCH --time=02:00:00
#SBATCH --nodes=2
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --partition=haswell
#SBATCH --reservation=READEX
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2500M
#SBATCH --comment="no_monitoring"
#SBATCH -J "PTFfoam"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ondrej.vysocky@vsb.cz

APP="simpleFoam -parallel"
THRDS=1
MPI_PROCS=28
PHASE_REG_NAME="iteration"

export FM_DIR=$( dirname "${BASH_SOURCE[0]}" )
if [ "$FM_DIR" == "." ]
then
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

export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_PLUGINS=hdeem_sync_plugin
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN="*/E"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_MPI_ENABLE_GROUPS=ENV

psc_frontend --apprun="$APP" --mpinumprocs=$MPI_PROCS --ompnumthreads=$THRDS --phase=$PHASE_REG_NAME --tune=readex_intraphase --config-file=readex_config_extended.xml --force-localhost --info=2 --selective-info=AutotuneAll,AutotunePlugins

