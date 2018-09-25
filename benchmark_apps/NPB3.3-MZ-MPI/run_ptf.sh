#!/bin/sh

#SBATCH --time=20:00:00   # walltime
#SBATCH --nodes=2  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
##SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH --mem=60000
#SBATCH -J "NPB-bt-mz_dta"   # job name
#SBATCH -A p_readex


#source ./env/set_env_npb_ptf.sh
source readex_env/set_env_ptf_hdeem.source
module load gdb

for i in `seq 0 23`
do
    cpufreq-set -c $i -f 2.5GHz
done

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin"
export SCOREP_RRL_VERBOSE="WARN"

#module load scorep-hdeem/sync-xmpi-gcc6.3
module list

export SCOREP_METRIC_PLUGINS=hdeem_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_MPI_ENABLE_GROUPS=ENV
export OMP_NUM_THREADS=12

# uncomment following line if using ATP library
# export ATP_EXECUTION_MODE=DTA
#gdbserver host:2345
PHASE=phase
app=bt-mz.C.2_ptf
psc_frontend --apprun="bin/$app" \
             --mpinumprocs=1 \
             --ompnumthreads=12 \
             --phase=$PHASE \
             --tune=readex_tuning \
             --config-file=ptf_readex_config.xml \
             --force-localhost \
             --info=7 \
             --selective-info=AutotuneAll,AutotunePlugins,FrontendStateMachines,AutotuneAgentStrategy \

