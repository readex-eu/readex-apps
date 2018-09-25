#!/bin/bash

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "bem4i_plain"   # job name
#SBATCH --output=bem4i_plain.out
#SBATCH --error=bem4i_plain.out
###############################################################################

module purge
source ./readex_env/set_env_plain.source
source ./set_env_bem4i.source

srun -N 1 -n 1 -c 24 ./dist/release_taurus/bem4i_plain #./input/cube_12.txt 5 1 4 4
 
