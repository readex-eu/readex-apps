#!/bin/sh

#SBATCH --time=5:00:00   # walltime
#SBATCH --nodes=2 # number of processor cores (i.e. tasks)
#SBATCH --ntasks=4
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -J "reg_rrl"   # job name
#SBATCH -A p_readex

source ../../init.sh
cd ${REL_PATH_APP_EXECUTION}

module purge
module use /projects/p_readex/modules
module load readex/ci_readex_bullxmpi1.2.8.4_gcc6.3.0

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_VERBOSE="VERBOSE"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin,OpenMPTP"
export SCOREP_RRL_TMM_PATH="tuning_model.json"
export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=false
export SCOREP_ENV_ENABLE_GROUPS=ENV

INPUT_FILE=in3.data

srun ./miniMD_openmpi_ptf -i ${INPUT_FILE}
