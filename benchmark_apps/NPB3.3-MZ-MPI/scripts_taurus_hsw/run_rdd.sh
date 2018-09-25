#!/bin/bash

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment=no_monitoring
#SBATCH --reservation=READEX
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J dyn_detect   # job name
#SBATCH --output=dyn_detect.out
#SBATCH --error=dyn_detect.out
###############################################################################
ps aux | grep diamon
echo "run RDD begin."

cd ..

module  purge
source readex_env/set_env_rdd.source

working_dir=$(pwd)
export LC_ALL=C

RDD_t=0.001 #1 ms
#RDD_t=0.01 #10 ms
#RDD_t=0.03 #30 ms
#RDD_t=0.1 #100 ms
#RDD_t=0.5 #500 ms
RDD_p=phase
RDD_c=0
RDD_v=0
RDD_w=0

export OMP_NUM_THREADS=12
export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
export SCOREP_FILTERING_FILE=scorep.filt 
export SCOREP_FILTERING_FILE="${working_dir}/scorep.filt"
export SCOREP_EXPERIMENT_DIRECTORY="${working_dir}/readex_dyn_detect"

app=bt-mz.C.2_rdd
rm -rf scorep-*
rm -f readex_config.xml

srun ./bin/$app

readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w "$SCOREP_EXPERIMENT_DIRECTORY/profile.cubex"

echo "RDD result = $?"

echo "run RDD done."
