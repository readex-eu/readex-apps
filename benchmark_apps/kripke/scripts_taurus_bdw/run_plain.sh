#!/bin/sh

#SBATCH --nodes=1
#SBATCH --tasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:30:00
#SBATCH -p broadwell
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH --comment="no_monitoring"

cd ../build
. ../readex_env/set_env_plain.source
. ../environment.sh

srun -n 28 ./kripke $KRIPKE_COMMAND

