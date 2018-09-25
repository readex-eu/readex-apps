#!/bin/sh

#SBATCH --nodes=1
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:30:00
#SBATCH -p haswell
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH --mem=60000
#SBATCH --comment="no_monitoring"

cd ../build
. ../readex_env/set_env_rrl.source
. ../environment.sh

#cp ../RESULTS/tuning_model.json .

export SCOREP_SUBSTRATE_PLUGINS='rrl'
export SCOREP_RRL_VERBOSE="WARN"
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin
export SCOREP_RRL_TMM_PATH=tuning_model.json
export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=false
export SCOREP_MPI_ENABLE_GROUPS=ENV

echo "running for RRL"
srun -n 24 ./kripke $KRIPKE_COMMAND

clearHdeem
startHdeem
srun -n 24 ./kripke $KRIPKE_COMMAND
stopHdeem
checkHdeem

echo "running for RRL done"
