cd ..

module purge
source ./readex_env/set_env_ptf_rapl.source

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

make veryclean
make -j 12

cp ./test/amg2013 ./test/amg2013_ptf

