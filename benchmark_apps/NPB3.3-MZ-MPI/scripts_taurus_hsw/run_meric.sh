#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=READEX
#SBATCH -J "NPB_meric"   # job name
#SBATCH --output=NPB_meric.out
##SBATCH --error=NPB_meric.out
###############################################################################

cd ..

# LOAD MODULES
module purge
source readex_env/set_env_meric.source
module list

export OMP_NUM_THREADS=12

# SET MERIC OUTPUT
#export MERIC_COUNTERS=perfevent
export MERIC_MODE=2
export MERIC_DEBUG=0
export MERIC_CONTINUAL=1
export MERIC_AGGREGATE=1
export MERIC_DETAILED=0

# FIRST RUN SOMETIMES CAUSES TROUBLE
export MERIC_FREQUENCY=25
export MERIC_UNCORE_FREQUENCY=30
export MERIC_NUM_THREADS=12
export MERIC_OUTPUT_DIR="bt-mz.C.2_meric_dir"

srun ./bin/bt-mz.C.2_meric
rm -rf $MERIC_OUTPUT_DIR
rm -rf ${MERIC_OUTPUT_DIR}Counters
#exit

echo $MERIC_OUTPUT_DIR

# FOR EACH SETTINGS
for thread in {12..4..-4}
do
  for cpu_freq in {25..13..-4}
  do
    for uncore_freq in {30..14..-4}
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
          srun --cpu_bind=verbose,sockets --nodes 1 --ntasks-per-node 2 --cpus-per-task 12 ./bin/bt-mz.C.2_meric
          ret=$?
        done
      done
    done
  done
done


~                                        
