#!/bin/sh

#source readex_env/set_env_ptf_hdeem.source
cd ..
source readex_env/set_env_ptf.source

#export SCOREP_PREP="scorep --online-access --user --nocompiler --thread=none --instrument-filter="${pathToFilter}"/scorep.filt"
#export SCOREP_PREP="scorep --online-access --user --nocompiler --nomemory --mpp=mpi --thread=none"
export SCOREP_PREP="scorep --online-access --user --nocompiler --nomemory --mpp=mpi --thread=omp"
export DSWITCH="-DUSE_SCOREP_MANUAL"

make clean
rm -rf lulesh2.0 lulesh2.0_ptf
make
mv lulesh2.0 lulesh2.0_ptf
