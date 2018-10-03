#!/bin/sh

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=6
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2500M
#SBATCH -J "READEX-espreso"
#SBATCH -A p_readex
	##SBATCH --mail-user=ondrej.vysocky@vsb.cz
#SBATCH --mail-type=ALL

################################################################################
source env/readex_env/set_env_plain.source
source env/modules.taurus

source env/paths.default
source env/threading.default 6

# first run without measurement
srun -N 1 -n 4 -c 6 --exclusive -p interactive --mem-per-cpu 2500M ./espreso_plain -c espreso_atp.ecf

# measurement
clearHdeem
startHdeem
srun -N 1 -n 4 -c 6 --exclusive -p interactive --mem-per-cpu 2500M ./espreso_plain -c espreso_atp.ecf
stopHdeem
checkHdeem


