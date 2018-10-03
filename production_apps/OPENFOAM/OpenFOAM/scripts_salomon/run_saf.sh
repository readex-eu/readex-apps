#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qexp
#PBS -N FOAM_SAF
#PBS -l select=1:ncpus=24:mpiprocs=24:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=01:00:00
#PBS -m be

APP='mpirun -n 24 simpleFoam -parallel'

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
export FM_DIR=$FM_DIR/..

cd $FM_DIR
source readex_env/set_env_saf.source
source scripts_$READEX_MACHINE/environment.sh

cd $FM_DIR/OpenFOAM-v1612+/
source $FM_DIR/OpenFOAM-v1612+/etc/bashrc
cd ../../motorBike24/

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

iteration=0
result=1
while [ $result != 0 ]; do
  iteration=$(($iteration +1))
  echo "result = "$result
  # run the application
  $APP
  echo "ITERATION done."
  $FM_DIR/scripts_$READEX_MACHINE/do_scorep_autofilter_single.sh 0.1
  result=$?
  echo "scorep_autofilter_singe done ($result)."
  mv scorep-* $iteration"_scorep"
  cp scorep.filt $iteration"_scorep.filt"
done
echo "end."
