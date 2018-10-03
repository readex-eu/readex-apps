#!/bin/sh

source ../../init.sh
cd ${REL_PATH_APP_EXECUTION}

#module purge
#module use /projects/p_readex/modules
#module load readex/ci_readex_bullxmpi1.2.8.4_gcc6.3.0

source env/modules.taurus.atp
source env/paths.default
source env/threading.default 12

INPUT_FILE=espreso.ecf

RDD_t=0.1
RDD_p=Main
RDD_c=10
RDD_v=10
RDD_w=10

export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm readex_config.xml

srun -N 1 -n 2 --ntasks-per-node=2 -c 12 --exclusive -p haswell --mem-per-cpu 2500M ./espreso_rdd

readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w scorep-*/profile.cubex

exit $?
