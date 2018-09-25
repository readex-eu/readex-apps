#!/bin/sh

#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2500M
#SBATCH --reservation=READEX
#SBATCH -A p_readex


LOC=$(pwd)
cd ..


module purge
source ./readex_env/set_env_rrl.source

export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin"
export SCOREP_RRL_TMM_PATH="tuning_model.json"
export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_RRL_CHECK_IF_RESET="no_reset"
export CHECK_IF_NODE_FULLY_OCCUPIED=0
export SCOREP_TUNING_CPU_FREQ_PLUGIN_VERBOSE=DEBUG
export SCOREP_TUNING_UNCORE_FREQ_PLUGIN_VERBOSE=DEBUG
clearHdeem
startHdeem
# run application
start_time=$(($(date +%s%N)/1000000))

srun  ./bin/bt-mz.B.24_ptf
stop_time=$(($(date +%s%N)/1000000))
stopHdeem
sleep 5
checkHdeem

time_total=$(echo "$stop_time - $start_time" | bc -l)
echo "Time=$time_total"

