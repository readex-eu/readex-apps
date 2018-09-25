#!/bin/sh

#SBATCH --time=01:00:00   # walltime
#SBATCH --nodes=4 
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem-per-cpu=2500M  
#SBATCH -A p_readex
#SBATCH -J "READEX-run_rdd"
#SBATCH --reservation=READEX

module load mkl/2017

SCRIPT_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
. ${SCRIPT_DIR}/readex_env/set_env_rdd.source
. ${SCRIPT_DIR}/../paths.source


export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
export SCOREP_FILTERING_FILE=scorep.filt

export NPROCS=96
export RDD_TIME=0.1

export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"

cd ${PROBLEM_PATH}

rm -f case.sif
rm -rf scorep-*

cp ${ELMER_ROOT}/case_cg.sif case.sif
cp ${ELMER_ROOT}/scorep.filt .

echo "Generating Grid for Elmer..."
if [ ! -d ${PROBLEM_PATH}/winkel/partitioning.96 ]; then
    echo "Problem winkel generated."
    ${ELMER_ROOT}/gen_problem.sh
else
    echo "Problem winkel already exists."
fi

echo "Running ElmerSolver for readex-dyn-detect..."
srun -n ${NPROCS} ${ELMER_HOME}/bin/ElmerSolver_mpi
echo "Running  ElmerSolver done."

echo "Running readex-dyn-detect..."
readex-dyn-detect -p "ElmerSolver_core" -t "${RDD_TIME}" -c 10 -v 10 -w 10 scorep-*/profile.cubex
#readex-dyn-detect -t "${RDD_TIME}" -c 10 -v 10 -w 10 scorep-*/profile.cubex

echo "RDD result = $?"

echo "Running readex-dyn-detect done." 

#cp -R scorep-* ../RESULTS/
cp readex_config.xml ${ELMER_ROOT}
