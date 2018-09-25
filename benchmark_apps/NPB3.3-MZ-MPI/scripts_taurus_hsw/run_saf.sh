#!/bin/sh

#SBATCH --time=6:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --reservation=READEX
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "NPB_C_autofilter"   # job name
#SBATCH --output=NPB_C_autofilter.out
###############################################################################

working_dir=$(pwd)
cd ..

module purge
source readex_env/set_env_saf.source
export LC_ALL=C

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_TOTAL_MEMORY=3G
export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM

export SCOREP_FILTERING_FILE="scorep.filt"
export SCOREP_EXPERIMENT_DIRECTORY="${working_dir}/readex_scorep_autofilt"

pathToFilter=$(pwd)

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

if [ "$READEX_INTEL" == "1" ]; then
  rm -rf old_scorep_icc.filt
  echo "" > scorep_icc.filt
fi

export OMP_NUM_THREADS=12

#SAF_t=0.001 #1 ms
SAF_t=0.1 #10 ms
#SAF_t=0.1 #100 ms
app=bt-mz.C.2_saf
result=1
while [ $result != 0 ]; do
  echo "result = "$result
  srun  ./bin/$app
  echo "Aplication run - done."
  ${working_dir}/do_scorep_autofilter_single.sh $SAF_t
  result=$?

if [ "$READEX_INTEL" == "1" ] && [ $result != 0 ]; then
  export FILTER_INTEL="-tcollect-filter="${pathToFilter}"/scorep_icc.filt"
  cd ${working_dir}
  ./compile_for_saf.sh
  cd ..
fi

  echo "scorep_autofilter_singe done ($result)."
done

echo "end of scorep-autofilter."

