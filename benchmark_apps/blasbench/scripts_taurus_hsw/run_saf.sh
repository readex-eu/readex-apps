#!/bin/bash

#SBATCH --time=1:00:00 # walltime
#SBATCH --nodes=1 # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2500M # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=READEX
#SBATCH -J "blasbench_saf" # job name
#SBATCH --output=blasbench_saf.out
#SBATCH --error=blasbench_saf.out
###############################################################################

LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_saf.source
source $LOC/set_env_blasbench.source

export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

if [ "$READEX_INTEL" == "1" ]; then
  rm -rf old_scorep_icc.filt
  echo "" > scorep_icc.filt
fi

#SAF_t=0.001 #1 ms
#SAF_t=0.01 #10 ms
SAF_t=0.1 #100 ms

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  srun --cpu_bind=verbose,sockets --nodes 1 --ntasks-per-node 2 --cpus-per-task 12 ./blasbench_saf --dgemv 300,8192 --dgemm 150,2048 --tan 50000 --phase 2
  echo "Aplication run - done."
  $LOC/do_scorep_autofilter_single.sh $SAF_t
  result=$?
 
if [ "$READEX_INTEL" == "1" ] && [ $result != 0 ]; then
  export FILTER_INTEL="-tcollect-filter=scorep_icc.filt"
  cd $LOC
  ./compile_for_saf.sh
  cd ..
fi

  echo "scorep_autofilter_single done ($result)."
done

echo "end of scorep-autofilter."
