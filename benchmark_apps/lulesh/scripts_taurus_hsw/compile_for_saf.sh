#!/bin/sh

cd ..
source readex_env/set_env_saf.source

#export SCOREP_PREP="scorep --online-access --user --compiler --thread=none"
export SCOREP_PREP="scorep --online-access --user --thread=omp --noopenmp --nocompiler --nomemory"
export DSWITCH="-DUSE_SCOREP"

make clean
rm -rf lulesh2.0 lulesh2.0_saf
make
mv lulesh2.0 lulesh2.0_saf
