#!/bin/sh

#SBATCH --nodes=4
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH -p haswell
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH --mem=60000
#SBATCH --reservation=READEX


module load mkl/2017
module load papi/5.5.1

SCRIPT_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
. ${SCRIPT_DIR}/readex_env/set_env_meric.source
. ${SCRIPT_DIR}/readex_env/set_env_plain.source
. ${SCRIPT_DIR}/../paths.source

cd ${ELMER_ROOT}

export NPROCS=96

export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"

# 2) Run
cd ${PROBLEM_PATH}
echo $PWD
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"
echo "SOLVER: ${SOLVER}"
rm -f case.sif
cp ${ELMER_ROOT}/case_cg.sif case.sif


stopHdeem
clearHdeem
startHdeem
sleep 1
stopHdeem


echo "Running Elmer (plain)"

clearHdeem
startHdeem

if [ ! -d ${PROBLEM_PATH}/winkel/partitioning.96 ]; then
    ${ELMER_ROOT}/gen_problem.sh
fi

srun -n $NPROCS ElmerSolver_mpi

stopHdeem
sleep 1
checkHdeem

echo "Running Elmer (plain) done"

