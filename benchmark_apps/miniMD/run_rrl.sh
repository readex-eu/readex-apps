#!/bin/sh

#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2500M
#SBATCH -J "miniMD_rrl"
#SBATCH -A p_readex

INPUT_FILE=in3.data #in.lj.miniMD
#####

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

# start RRL-tuned run
module purge
source ./readex_env/set_env_rrl.source

export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin,OpenMPTP"
export SCOREP_RRL_TMM_PATH="tuning_model.json"
export SCOREP_MPI_ENABLE_GROUPS=ENV

# run RRL-tuned application
srun ./miniMD_openmpi_ptf -i $INPUT_FILE
# end RRL-tuned run
