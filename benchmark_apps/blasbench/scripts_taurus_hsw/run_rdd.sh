#!/bin/bash

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=READEX
#SBATCH -J "blasbench_rdd"   # job name
#SBATCH --output=blasbench_rdd.out
#SBATCH --error=blasbench_rdd.out
###############################################################################

echo "run RDD begin."

LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_rdd.source
source $LOC/set_env_blasbench.source

#RDD_t=0.001 #1 ms
#RDD_t=0.01 #10 ms
#RDD_t=0.03 #30 ms
RDD_t=0.1 #100 ms
#RDD_t=0.5 #500 ms
RDD_p=Main
RDD_c=10
RDD_v=10
RDD_w=10

export SCOREP_PROFILING_FORMAT=cube_tuple
#export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM

rm -rf scorep-*
rm -f readex_config.xml

srun --cpu_bind=verbose,sockets --nodes 1 --ntasks-per-node 2 --cpus-per-task 12 ./blasbench_rdd --dgemv 300,8192 --dgemm 150,2048 --tan 60000 --phase 2

readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w scorep-*/profile.cubex

echo "RDD result = $?"

echo "run RDD done."
