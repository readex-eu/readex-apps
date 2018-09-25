#!/bin/bash

#SBATCH --time=0:30:00   # walltime
#SBATCH --nodes=4  # number of processor cores (i.e. tasks)
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=14
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2200M   # memory per CPU core
#SBATCH -A p_readex
##SBATCH --reservation=p_readex_56
#SBATCH -J "amg2013_plain"   # job name
#SBATCH --output=amg2013_plain.out
#SBATCH --error=amg2013_plain.out
###############################################################################

cd ..

module purge
source ./readex_env/set_env_plain.source

srun --cpu_bind=verbose,sockets --nodes 4 --ntasks-per-node 2 --cpus-per-task 14 ./test/amg2013_plain -P 2 2 2 -r 40 40 40 

