#!/bin/bash

#SBATCH --time=0:30:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=14
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --mem-per-cpu=2000M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=p_readex_56
##SBATCH -J "NPB_bt_plain"   # job name
##SBATCH --output=NPB_bt_plain.out
##SBATCH --error=NPB_bt_plain.out

###############################################################################

cd ..
source readex_env/set_env_plain.source
module list

export OMP_NUM_THREADS=14

clearHdeem
startHdeem
# run application
start_time=$(($(date +%s%N)/1000000))
app=bt-mz.C.2_plain
cd bin
srun $app
stop_time=$(($(date +%s%N)/1000000))
stopHdeem
sleep 5
checkHdeem

time_total=$(echo "$stop_time - $start_time" | bc -l)
echo "Time=$time_total"

