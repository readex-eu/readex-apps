#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "bem4i_meric"   # job name
#SBATCH --output=bem4i_meric.out
#SBATCH --error=bem4i_meric.out
###############################################################################

# LOAD MODULES
module purge
source ./readex_env/set_env_meric.source
source ./set_env_bem4i.source

# SET MERIC OUTPUT
#export MERIC_COUNTERS=perfevent
export MERIC_MODE=0
export MERIC_DEBUG=0
export MERIC_CONTINUAL=1
export MERIC_AGGREGATE=0
export MERIC_DETAILED=0

#NUMA=numactl" "--interleave=all"
NUMA=
PRG=srun" "-N" "1" "-n" "1" "-c" "24" "$NUMA" "./dist/release_taurus/bem4i_meric

# FIRST RUN SOMETIMES CAUSES TROUBLE
export MERIC_NUM_THREADS=24
export MERIC_FREQUENCY=25
export MERIC_UNCORE_FREQUENCY=30
export MERIC_OUTPUT_DIR="bem4i_meric_dir"
$PRG
rm -r $MERIC_OUTPUT_DIR
rm -r ${MERIC_OUTPUT_DIR}Counters
#exit

echo $MERIC_OUTPUT_DIR

# FOR EACH SETTINGS
for thread in {24..4..-4}
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
          $PRG              
          ret=$?
        done
      done 
    done
  done
done

