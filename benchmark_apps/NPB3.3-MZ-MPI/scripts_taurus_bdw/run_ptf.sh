#!/bin/bash

#SBATCH --time=08:00:00   # walltime
#SBATCH --nodes=2  # number of processor cores (i.e. tasks)
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=14
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --reservation=p_readex_56
#SBATCH --mem-per-cpu=2000M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "bt_C_dta"   # job name
#SBATCH --output=bt_C_dta.out
#SBATCH --error=bt_C_dta.out
###############################################################################

cd ..

module purge
source readex_env/set_env_ptf_rapl.source
#module load scorep-uncore
#module load scorep-uncore/sync/2016-09-13
#module load scorep-apapi
#module load intel
module list

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

echo $LD_LIBRARY_PATH

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin,OpenMPTP
export SCOREP_RRL_VERBOSE="INFO"
#export SCOREP_RRL_SIGNIFICANT_DURATION_MS=10
#export SCOREP_TOTAL_MEMORY=3G
export SCOREP_TUNING_CPU_FREQ_PLUGIN_VERBOSE=DEBUG
export SCOREP_TUNING_UNCORE_FREQ_PLUGIN_VERBOSE=DEBUG
export SCOREP_METRIC_PLUGINS="x86_energy_sync_plugin"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN="*/E"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE=DEBUG
export SCOREP_METRIC_PLUGINS_SEP=";"

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true


#export SCOREP_METRIC_PLUGINS="upe_plugin,apapi_plugin"
#export SCOREP_METRIC_SCOREP_SUBSTRATE_RRL="cpu_freq_plugin,uncore_freq_plugin"
#export SCOREP_METRIC_UPE_PLUGIN="hswep_unc_cbo0::UNC_C_CLOCKTICKS"
#export SCOREP_METRIC_PAPI="PAPI_TOT_CYC"


export OMP_NUM_THREADS=14

PHASE=phase
app=bt-mz.C.2_ptf
psc_frontend --apprun="./bin/bt-mz.C.2_ptf" --mpinumprocs=2 --ompnumthreads=14 --phase=$PHASE --tune=readex_intraphase --config-file=readex_config.xml --force-localhost --info=2 --selective-info=AutotuneAll,AutotunePlugins

