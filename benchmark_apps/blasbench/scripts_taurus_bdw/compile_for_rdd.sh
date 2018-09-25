LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_rdd.source
source $LOC/set_env_blasbench.source

if [ "$READEX_INTEL" == "1" ]; then
  export FILTER_GNU=""
  export FILTER_INTEL="-tcollect-filter=scorep_icc.filt"
else
  export FILTER_GNU="--instrument-filter=scorep.filt"
  export FILTER_INTEL=""
fi

make clean && make dyndetect
cp blasbench blasbench_rdd

