#!/bin/sh

#SBATCH --nodes=1
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:30:00
#SBATCH -p haswell
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH --mem=60000
#SBATCH --comment="no_monitoring"

cd ../build
. ../readex_env/set_env_plain.source
. ../environment.sh

stopHdeem
clearHdeem
startHdeem
sleep 1
stopHdeem


echo "running kripke"

clearHdeem
startHdeem
srun -n 24 ./kripke $KRIPKE_COMMAND
stopHdeem
sleep 1
checkHdeem

echo "running kripke done"
