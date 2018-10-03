#!/bin/sh

#SBATCH --time=01:00:00
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


#####################################################################
cd ../../

source env/readex_env/set_env_plain.source
source env/environment.sh

source env/paths.default
source env/threading.default 14

# first run without measurement
srun -N 1 -n 2 -c 14 --exclusive -p interactive --mem-per-cpu 2200M ./espreso_plain -c espreso.ecf



