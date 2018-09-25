#!/bin/bash

#SBATCH --time=0:30:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=READEX
#SBATCH -J "blasbench_plain"   # job name
#SBATCH --output=blasbench_plain.out
#SBATCH --error=blasbench_plain.out
###############################################################################

LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_plain.source
source $LOC/set_env_blasbench.source

srun --cpu_bind=verbose,sockets --nodes 1 --ntasks-per-node 2 --cpus-per-task 12 ./blasbench_plain --dgemv 300,8192 --dgemm 150,2048 --tan 60000 --phase 10 

