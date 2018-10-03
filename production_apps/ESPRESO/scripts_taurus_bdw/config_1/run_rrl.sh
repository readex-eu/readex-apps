#!/bin/sh

#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ondrej.vysocky@vsb.cz

#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=14

#SBATCH --partition=broadwell 
#SBATCH --reservation=p_readex_44
#SBATCH -A p_readex	#to account your compute time on the readex project
#SBATCH --exclusive	
#SBATCH --mem-per-cpu=2200M
#SBATCH --comment="no_monitoring"

cd ../../
source env/readex_env/set_env_rrl.source
source env/environment.sh

source env/paths.default
source env/threading.default 14

export SCOREP_SUBSTRATE_PLUGINS='rrl'
export SCOREP_RRL_VERBOSE="WARN"
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin,OpenMPTP
export SCOREP_RRL_TMM_PATH=tuning_model.json
export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=false
export SCOREP_MPI_ENABLE_GROUPS=ENV

echo "running for RRL"
srun -n 2 -c 14 ./espreso_rdd -c espreso.ecf


echo "running for RRL done"
