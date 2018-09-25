#!/bin/sh

cd ..
module purge
source readex_env/set_env_rdd.source


#export SCOREP_PREP="scorep --online-access --user --nocompiler --thread=none --instrument-filter="${pathToFilter}"/scorep.filt"
export SCOREP_PREP="scorep --online-access --user --nocompiler --thread=omp --noopenmp"
export DSWITCH="-DUSE_SCOREP_MANUAL"

make clean
rm -rf lulesh2.0 lulesh2.0_rdd
make
mv lulesh2.0 lulesh2.0_rdd
