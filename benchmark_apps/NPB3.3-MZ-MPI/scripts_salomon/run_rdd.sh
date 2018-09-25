#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N BTMZ_rdd
#PBS -o BTMZ_rdd.out
#PBS -e BTMZ_rdd.err
#PBS -l select=1:ncpus=24:mpiprocs=2:ompthreads=12:accelerator=true
#PBS -l x86_adapt=true
#PBS -l walltime=01:00:00
###############################################################################

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
  LOC=$PBS_O_WORKDIR
else
  LOC=$(pwd)
fi

cd $LOC/..

TMPFILE=btmz_rdd.tmp
rm -f $TMPFILE

module  purge
source readex_env/set_env_rdd.source

working_dir=$(pwd)
export LC_ALL=C

RDD_t=0.0001 #1 ms
#RDD_t=0.01 #10 ms
#RDD_t=0.03 #30 ms
#RDD_t=0.1 #100 ms
#RDD_t=0.5 #500 ms
RDD_p=phase
RDD_c=0
RDD_v=0
RDD_w=0

export SCOREP_PROFILING_FORMAT=cube_tuple
#export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
export SCOREP_EXPERIMENT_DIRECTORY="readex_dyn_detect"

app=bt-mz.C.2_rdd
rm -rf readex_dyn_detect-*
rm -f readex_config.xml

mpirun $BIND_TO_SOCKETS ./bin/$app

readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w "$SCOREP_EXPERIMENT_DIRECTORY/profile.cubex"

echo "RDD result = $?"

echo "run RDD done."
