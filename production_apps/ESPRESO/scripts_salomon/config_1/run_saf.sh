#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qexp
#PBS -N espreso_saf
#PBS -l select=1:ncpus=24:mpiprocs=2:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=01:00:00  
#PBS -m be

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
cd $FM_DIR/../../

source env/readex_env/set_env_saf.source
source env/environment.sh

source env/paths.default
source env/threading.default 12


export SCOREP_TOTAL_MEMORY=3G
export SCOREP_FILTERING_FILE=scorep.filt

#SAF_t=0.001 #1 ms
SAF_t=0.01 #10 ms
#SAF_t=0.1 #100 ms

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

if [ "$READEX_INTEL" == "1" ]; then
  rm -rf old_scorep_icc.filt
  echo "" > scorep_icc.filt
fi

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  mpirun -n 2 ./espreso_saf -c espreso_scaling.ecf
  echo "Aplication run - done."
  ./do_scorep_autofilter_single.sh $SAF_t
  result=$?
  echo "scorep_autofilter_singe done ($result)."
done

if [ "$READEX_INTEL" == "1" ]
then
	echo ".*:*__gnu_cxx::* off" >> scorep_icc.filt
	echo ".*:*std::* off" >> scorep_icc.filt
fi
echo "end of scorep-autofilter."
