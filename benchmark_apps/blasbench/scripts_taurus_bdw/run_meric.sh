#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=14
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2200M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=p_readex_56
#SBATCH -J "blasbench_meric"   # job name
#SBATCH --output=blasbench_meric.out
#SBATCH --error=blasbench_meric.out
###############################################################################

LOC=$(pwd)
cd ..

# LOAD MODULES
module purge
source ./readex_env/set_env_meric.source
source $LOC/set_env_blasbench.source

# SET MERIC OUTPUT
#export MERIC_COUNTERS=perfevent
export MERIC_MODE=1
export MERIC_DEBUG=0
export MERIC_CONTINUAL=1
export MERIC_AGGREGATE=1
export MERIC_DETAILED=0

PRG=srun" "--cpu_bind=verbose,sockets" "--nodes" "1" "--ntasks-per-node" "2" "--cpus-per-task" "14" "./blasbench_meric" "--dgemv" "300,8192" "--dgemm" "150,2048" "--tan" "60000" "--phase" "2

# FIRST RUN SOMETIMES CAUSES TROUBLE
export MERIC_NUM_THREADS=14
export MERIC_FREQUENCY=24
export MERIC_UNCORE_FREQUENCY=27
export MERIC_OUTPUT_DIR="blasbench_meric_dir"
$PRG
rm -rf $MERIC_OUTPUT_DIR
rm -rf ${MERIC_OUTPUT_DIR}Counters
#exit

echo $MERIC_OUTPUT_DIR

# FOR EACH SETTINGS
for thread in {14..6..-4}
do
  for cpu_freq in {24..12..-4}
  do
    for uncore_freq in {27..15..-4}
    do
      export MERIC_OUTPUT_FILENAME=$thread"_"$cpu_freq"_"$uncore_freq
      export MERIC_NUM_THREADS=$thread
      export MERIC_FREQUENCY=$cpu_freq
      export MERIC_UNCORE_FREQUENCY=$uncore_freq

      echo
      echo TEST $MERIC_OUTPUT_FILENAME
      for repeat in {1..1..1}
      do
        ret=1
        while [ "$ret" -ne 0 ]
        do
          $PRG              
          ret=$?
        done
      done 
    done
  done
done

