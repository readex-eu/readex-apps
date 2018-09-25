#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N BTMZ_meric
#PBS -o BTMZ_meric.out
#PBS -e BTMZ_meric.err
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

TMPFILE=btmz_meric.tmp
rm -f $TMPFILE

# LOAD MODULES
module purge
source readex_env/set_env_meric.source
module list

export OMP_NUM_THREADS=12

# SET MERIC OUTPUT
#export MERIC_COUNTERS=perfevent
export MERIC_MODE=1
export MERIC_DEBUG=0
export MERIC_CONTINUAL=1
export MERIC_AGGREGATE=1
export MERIC_DETAILED=0

# FIRST RUN SOMETIMES CAUSES TROUBLE
export MERIC_FREQUENCY=25
export MERIC_UNCORE_FREQUENCY=30
export MERIC_NUM_THREADS=12
export MERIC_OUTPUT_DIR="bt-mz.C.2_meric_dir"

mpirun $BIND_TO_SOCKETS ./bin/bt-mz.C.2_meric
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
          mpirun $BIND_TO_SOCKETS ./bin/bt-mz.C.2_meric
          ret=$?
        done
      done
    done
  done
done

rm -f $TMPFILE
~                                        
