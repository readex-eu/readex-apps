#!/bin/sh

#SBATCH --time=20:00:00   # walltime
#SBATCH --nodes=2  # number of processor cores (i.e. tasks)
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=6
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2500M
#SBATCH -J "ATP-espreso"
#SBATCH -A p_readex
	##SBATCH --mail-user=ondrej.vysocky@vsb.cz
#SBATCH --mail-type=ALL


handle_error() {
        if [ $1 != 0 ]; then
                echo "Error from job_build.sh on $2."
                exit 1
        fi
}

source env/readex_env/set_env_ptf_hdeem.source
source env/modules.taurus

source env/paths.default
source env/threading.default 6

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin,OpenMPTP
export SCOREP_RRL_VERBOSE="WARN"

export SCOREP_METRIC_PLUGINS=hdeem_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN="*/E"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_MPI_ENABLE_GROUPS=ENV

export ATP_EXECUTION_MODE=DTA

PHASE=Main

psc_frontend --apprun="./espreso_atp -c espreso_atp.ecf" --mpinumprocs=4 --ompnumthreads=6 --phase=$PHASE --tune=readex_intraphase --config-file=readex_config.xml --force-localhost --info=1 --selective-info=AutotuneAll,AutotunePlugins,FrontendStateMachines,AutotuneAgentStrategy,ApplTuningParameter


