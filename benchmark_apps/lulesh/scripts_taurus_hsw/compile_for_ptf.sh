#!/bin/sh

source readex_env/set_env_ptf_hdeem.source
module avail readex

pathToFilter=$(pwd)

#export SCOREP_PREP="scorep --online-access --user --compiler --thread=none --instrument-filter="${pathToFilter}"/scorep.filt"

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
rm -rf lulesh2.0 lulesh2.0_ptf
make
mv lulesh2.0 lulesh2.0_ptf
