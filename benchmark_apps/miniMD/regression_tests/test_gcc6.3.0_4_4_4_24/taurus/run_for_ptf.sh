#!/bin/sh

#SBATCH --time=5:00:00   # walltime
#SBATCH --nodes=2 # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -J "reg_ptf"   # job name
#SBATCH -A p_readex

source ../../init.sh
cd ${REL_PATH_APP_EXECUTION}

module purge
module use /projects/p_readex/modules
module load readex/ci_readex_bullxmpi1.2.8.4_gcc6.3.0

INPUT_FILE=in3.data
PHASE=INTEGRATE_RUN_LOOP
NP=2 # check against --ntasks and tasks-per-node

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export LD_LIBRARY_PATH+=:/projects/p_readex/readex-atp/atp-2017-09-07-gcc6.3.0-bullxmpi/lib

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin,OpenMPTP
export SCOREP_RRL_VERBOSE="WARN"

#module load scorep-hdeem/sync-2017-01-31-git-hdeem2.2.20ms-xmpi-gcc5.3
module load scorep-hdeem/sync-2017-04-13-git-hdeem2.2.20ms-xmpi-gcc6.3
export SCOREP_METRIC_PLUGINS=hdeem_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_MPI_ENABLE_GROUPS=ENV

psc_frontend --apprun="./miniMD_openmpi_ptf -i $INPUT_FILE" --mpinumprocs=$NP --ompnumthreads=12 --phase=$PHASE --tune=readex_intraphase --config-file=ptf_readex_config.xml --force-localhost --info=7 #--selective-info=AutotuneAll,AutotunePlugins

