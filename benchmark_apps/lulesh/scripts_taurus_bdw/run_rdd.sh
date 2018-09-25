#!/bin/bash


#SBATCH --time=02:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mohak.chadha@tum.de

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=28

#SBATCH --partition=broadwell
#SBATCH --reservation=p_readex_56
#SBATCH -A p_readex 
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2000M   # memory per CPU core
#SBATCH --comment="no_monitoring"
#SBATCH --output=lulesh_rdd.out
#SBATCH --error=lulesh_rdd.out


cd ..

module purge
source ./readex_env/set_env_rdd.source

RDD_t=0.1
RDD_p=foo
RDD_c=10
RDD_v=10
RDD_w=10

export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
#export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm readex_config.xml

#export GOMP_CPU_AFFINITY=0-27
#srun -N 1 -n 1 -c 24 --exclusive -p haswell --mem-per-cpu 2500M ./lulesh2.0_rdd -i 50 -s 37
srun -n 1 -c 28 ./lulesh2.0_rdd -i 100 -s 150 
readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w scorep-*/profile.cubex

echo "RDD result = $?"

echo "run RDD done."

