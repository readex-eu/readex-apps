#!/bin/bash

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=4  # number of processor cores (i.e. tasks)
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=14
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2200M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH --reservation=p_readex_56
#SBATCH -J "amg2013_rdd"   # job name
#SBATCH --output=amg2013_rdd.out
#SBATCH --error=amg2013_rdd.out
###############################################################################

echo "run RDD begin."

cd ..

module purge
source ./readex_env/set_env_rdd.source

#RDD_t=0.001 #1 ms
#RDD_t=0.01 #10 ms
#RDD_t=0.03 #30 ms
RDD_t=0.1 #100 ms
#RDD_t=0.5 #500 ms
RDD_p=main_phase
RDD_c=10
RDD_v=10
RDD_w=10

export SCOREP_PROFILING_FORMAT=cube_tuple
#export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM

rm -rf scorep-*
rm -f readex_config.xml

srun --cpu_bind=verbose,sockets --nodes 4 --ntasks-per-node 2 --cpus-per-task 14 ./test/amg2013_rdd -P 2 2 2 -r 40 40 40

readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w scorep-*/profile.cubex

echo "RDD result = $?"

echo "run RDD done."
