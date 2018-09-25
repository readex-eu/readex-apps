#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N lulesh_saf
#PBS -o lulesh_saf.out
#PBS -e lulesh_saf.err
#PBS -l select=1:ncpus=24:ompthreads=24:accelerator=true
#PBS -l walltime=01:00:00
###############################################################################

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
  LOC=$PBS_O_WORKDIR
else
  LOC=$(pwd)
fi

cd $LOC/..

TMPFILE=lulesh_saf.tmp
rm -f $TMPFILE

module purge
source ./readex_env/set_env_saf.source

export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

if [ "$READEX_INTEL" == "1" ]; then
  echo ".*:*__gnu_cxx::* off" >> scorep_icc.filt
  echo ".*:*std::* off" >> scorep_icc.filt
  export FILTER_INTEL="-tcollect-filter=scorep_icc.filt"
  cd $LOC
  ./compile_for_saf.sh 2>&1 | tee -a ../$TMPFILE
  cd $LOC/..
fi

#SAF_t=0.001 #1 ms
#SAF_t=0.01 #10 ms
SAF_t=0.1 #100 ms

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  ./lulesh2.0_saf -i 100 -s 150 2>&1 | tee -a $TMPFILE
  echo "Aplication run - done."
  $LOC/do_scorep_autofilter_single.sh $SAF_t
  result=$?

if [ "$READEX_INTEL" == "1" ] && [ $result != 0 ]; then
  export FILTER_INTEL="-tcollect-filter=scorep_icc.filt"
  cd $LOC
  ./compile_for_saf.sh 2>&1 | tee -a ../$TMPFILE
  cd $LOC/..
fi

  echo "scorep_autofilter_single done ($result)."
done

rm -f $TMPFILE

