#!/bin/bash

#SBATCH --time=0:30:00   # walltime
#SBATCH --nodes=4  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "miniMD_plain"   # job name
#SBATCH --output=miniMD_plain.out
#SBATCH --error=miniMD_plain.err
###############################################################################

module purge
source /readex_env/set_env_plain.source

INPUT_FILE=in3.data

srun --cpu_bind=verbose,sockets --nodes 4 --ntasks-per-node 2 --cpus-per-task 12 ./miniMD_openmpi_plain -i $INPUT_FILE
