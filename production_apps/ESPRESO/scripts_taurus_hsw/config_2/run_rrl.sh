#!/bin/sh

#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=6
#SBATCH --time=0-00:30:00
#SBATCH -p haswell
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH -J "RRL-espreso"
#SBATCH --mem-per-cpu=2500M
#SBATCH --comment="no_monitoring"

source env/readex_env/set_env_ptf.source
source env/modules.taurus

source env/paths.default
source env/threading.default 6

export SCOREP_SUBSTRATE_PLUGINS='rrl'
export SCOREP_RRL_VERBOSE="WARN"
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin,OpenMPTP
export SCOREP_RRL_TMM_PATH=tuning_model.json
export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=false
export SCOREP_MPI_ENABLE_GROUPS=ENV

echo "running for RRL"
srun -n 4 -c 6 ./espreso_ptf -c espreso_atp.ecf

clearHdeem
startHdeem
srun -n 4 -c 6 ./espreso_ptf -c espreso_atp.ecf
stopHdeem
checkHdeem

echo "running for RRL done"
