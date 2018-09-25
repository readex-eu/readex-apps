#!/bin/bash

##SBATCH --time=1:00:00   # walltime
##SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=24
##SBATCH --exclusive
##SBATCH --partition=haswell
##SBATCH --comment="cpufreqchown"
##SBATCH --mem-per-cpu=2500M   # memory per CPU core
##SBATCH -A p_readex
##SBATCH -J $app   # job name
##SBATCH --output=$app.out
##SBATCH --error=$app.out
###############################################################################

echo "run RDD begin."

#source env/set_env_npb_ptf.sh
source readex_env/set_env_rdd.source

#RDD_t=0.001 #1 ms
#RDD_t=0.01 #10 ms
#RDD_t=0.03 #30 ms
RDD_t=0.001 #100 ms
#RDD_t=0.5 #500 ms
RDD_p=phase
RDD_c=10
RDD_v=10
RDD_w=10

export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
export SCOREP_FILTERING_FILE=scorep.filt 
app=bt-mz.C.2_rdd
rm -rf scorep-*
rm -f readex_config.xml

srun -N 1 -n 2 -c 12 --exclusive -p interactive --mem-per-cpu 2500M bin/$app

readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w scorep-*/profile.cubex

echo "RDD result = $?"

echo "run RDD done."
