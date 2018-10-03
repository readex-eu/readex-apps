#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N espreso_rdd
#PBS -l select=1:ncpus=24:mpiprocs=2:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=24:00:00  
#PBS -m be

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
cd $FM_DIR/../../

echo "run RDD begin."

source env/readex_env/set_env_rdd.source
source env/environment.sh

source env/paths.default
source env/threading.default 12

#RDD_t=0.001 #1 ms
#RDD_t=0.01 #10 ms
#RDD_t=0.03 #30 ms
RDD_t=0.1 #100 ms
#RDD_t=0.5 #500 ms
RDD_p=Main
RDD_c=10
RDD_v=10
RDD_w=10

export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM

rm -rf scorep-*
rm readex_config.xml

mpirun -n 2 ./espreso_rdd -c espreso.ecf

readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w scorep-*/profile.cubex

echo "RDD result = $?"

echo "run RDD done."
