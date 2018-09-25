#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N blasbench_meric
#PBS -o blasbench_meric.out
#PBS -e blasbench_meric.err
#PBS -l select=1:ncpus=24:mpiprocs=2:ompthreads=12:accelerator=true
#PBS -l x86_adapt=true
#PBS -l walltime=06:00:00
###############################################################################

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
  LOC=$PBS_O_WORKDIR
else
  LOC=$(pwd)
fi

cd $LOC/..

TMPFILE=blasbench_meric.tmp
rm -f $TMPFILE

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

# FIRST RUN SOMETIMES CAUSES TROUBLE
export MERIC_NUM_THREADS=12
export MERIC_FREQUENCY=25
export MERIC_UNCORE_FREQUENCY=30
export MERIC_OUTPUT_DIR=blasbench_meric_dir
mpirun $BIND_TO_SOCKETS ./blasbench_meric --dgemv 300,8192 --dgemm 150,2048 --tan 60000 --phase 2 
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
      echo TEST $MERIC_OUTPUT_FILENAME | tee -a $TMPFILE
      for repeat in {1..1..1}
      do
        ret=1
        while [ "$ret" -ne 0 ]
        do
          mpirun $BIND_TO_SOCKETS ./blasbench_meric --dgemv 300,8192 --dgemm 150,2048 --tan 60000 --phase 2 | tee -a $TMPFILE 
          ret=$?
        done
      done 
    done
  done
done

rm -f $TMPFILE

