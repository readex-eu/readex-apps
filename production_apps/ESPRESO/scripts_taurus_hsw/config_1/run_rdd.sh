#!/bin/sh

#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -J "RDD-espreso"   # job name
#SBATCH -A p_readex
	##SBATCH --mail-user=ondrej.vysocky@vsb.cz
#SBATCH --mail-type=ALL
#SBATCH --reservation=READEX

echo "run RDD begin."
cd ../../
source env/readex_env/set_env_rdd.source
source env/environment.sh

source env/paths.default
source env/threading.default 12

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
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM

rm -rf scorep-*
rm readex_config.xml

srun -n 2 -c 12 ./espreso_rdd -c espreso.ecf

readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w scorep-*/profile.cubex

echo "RDD result = $?"

echo "run RDD done."
