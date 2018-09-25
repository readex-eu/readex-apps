#!/bin/sh

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "amg2013_rdd"   # job name
#SBATCH --output=amg2013_rdd.out
#SBATCH --error=amg2013_rdd.out
###############################################################################

module purge
source ./readex_env/set_env_rdd.source

INPUT_FILE=in3.data #in.lj.miniMD

RDD_t=0.001
RDD_p=INTEGRATE_RUN_LOOP
RDD_c=10
RDD_v=10
RDD_w=10

export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm readex_config.xml

srun ./miniMD_openmpi_rdd -i $INPUT_FILE

readex-dyn-detect -t $RDD_t -p $RDD_p -c $RDD_c -v $RDD_v -w $RDD_w scorep-*/profile.cubex

echo "RDD result = $?"

echo "run RDD done."
