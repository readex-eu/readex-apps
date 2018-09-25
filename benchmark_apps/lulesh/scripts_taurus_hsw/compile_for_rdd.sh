#!/bin/sh

cd ..
source readex_env/set_env_rdd.source

pathToFilter=$(pwd)

#export SCOREP_PREP="scorep --online-access --user --compiler --thread=none --instrument-filter="${pathToFilter}"/scorep.filt"
#export SCOREP_PREP="scorep --online-access --user --compiler --mpp=mpi --thread=omp --noopenmp"
#export DSWITCH="-DUSE_SCOREP"


pathToFilter=$(pwd)

if [ "$READEX_INTEL" == "1" ]; then
	export FILTER_GNU=""
	export FILTER_INTEL="-tcollect-filter="${pathToFilter}"/scorep_icc.filt"
else
	export FILTER_GNU="--instrument-filter="${pathToFilter}"/scorep.filt"
	export FILTER_INTEL=""
fi

export SCOREP_PREP="scorep --online-access --user --mpp=mpi --thread=omp --noopenmp "${FILTER_GNU}
export DSWITCH="-DUSE_SCOREP"

make clean
rm -rf lulesh2.0 lulesh2.0_rdd
make
mv lulesh2.0 lulesh2.0_rdd
