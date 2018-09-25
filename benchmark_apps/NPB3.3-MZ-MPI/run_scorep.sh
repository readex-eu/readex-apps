#!/bin/sh

##SBATCH --time=1:00:00   # walltime
##SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
##SBATCH --ntasks-per-node=2
##SBATCH --cpus-per-task=12
##SBATCH --exclusive
##SBATCH --partition=haswell
##SBATCH --comment="cpufreqchown"
##SBATCH --mem-per-cpu=2500M   # memory per CPU core
##SBATCH -A p_readex
##SBATCH -J ""$app"_saf"   # job name
##SBATCH --output="$app"_saf.out
##SBATCH --error="$app"_saf.out
###############################################################################



source readex_env/set_env_saf.source

export SCOREP_FILTERING_FILE=scorep.filt

srun bin/bt-mz.C.2_saf

