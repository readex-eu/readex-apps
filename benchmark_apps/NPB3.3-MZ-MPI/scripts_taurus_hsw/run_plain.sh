#!/bin/bash

#SBATCH --time=02:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --reservation=READEX
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "NPB_plain"   # job name
#SBATCH --output=NPB_plain.out
###############################################################################

ps aux | grep diamon


cd ..

source readex_env/set_env_plain.source
module list

export OMP_NUM_THREADS=12

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


