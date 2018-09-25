#!/bin/bash
#SBATCH --time=02:30:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -J "lulesh"   # job name
#SBATCH -A p_readex
#SBATCH --output=lulesh_rdd.out
#SBATCH --error=lulesh_rdd.out 
#SBATCH --reservation=READEX

cd ..

module purge
source ./readex_env/set_env_rdd.source

RDD_t=0.1
RDD_p=foo
RDD_c=10
RDD_v=10
RDD_w=10

export SCOREP_PROFILING_FORMAT=cube_tuple
#export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
#export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm readex_config.xml

srun -N 1 -n 1 -c 24 --exclusive -p haswell --mem-per-cpu 2500M ./lulesh2.0_rdd -p -i 100 -s 150
readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w scorep-*/profile.cubex

echo "RDD result = $?"

echo "run RDD done."

