#!/bin/bash

#SBATCH --time=30:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=14
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --mem-per-cpu=2000M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=p_readex_56
#SBATCH -J "bt_meric"   # job name
#SBATCH --output=bt_meric.out
#SBATCH --error=bt_meric.out
###############################################################################

cd ..
# LOAD MODULES
module purge
source readex_env/set_env_meric.source
module list

# SET MERIC OUTPUT
#export MERIC_COUNTERS=perfevent
export MERIC_MODE=2
export MERIC_DEBUG=0
export MERIC_CONTINUAL=1
export MERIC_AGGREGATE=1
export MERIC_DETAILED=0

# FIRST RUN SOMETIMES CAUSES TROUBLE
export MERIC_NUM_THREADS=14
export MERIC_FREQUENCY=24
export MERIC_UNCORE_FREQUENCY=27
export MERIC_OUTPUT_DIR="bt-mz.C.2_meric_dir"
srun ./bin/bt-mz.C.2_meric
rm -rf $MERIC_OUTPUT_DIR
rm -rf ${MERIC_OUTPUT_DIR}Counters
#exit

echo $MERIC_OUTPUT_DIR

# FOR EACH SETTINGS
for thread in {14..4..-4}
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
          srun ./bin/bt-mz.C.2_meric
          ret=$?
        done
      done
    done
  done
done

          
