#!/bin/sh

#SBATCH -t 01:0:00
#SBATCH --mem=24000
#SBATCH --nodes=4
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH -p haswell
#SBATCH --comment="no_monitoring"
#SBATCH --reservation=READEX


# 1) Modules and variables
module load mkl/2017
module load papi/5.5.1

SCRIPT_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
. ${SCRIPT_DIR}/readex_env/set_env_meric.source
. ${SCRIPT_DIR}/../paths.source

cd ${ELMER_ROOT}


export LD_LIBRARY_PATH+=:/projects/p_readex/it4i/mericGCC/lib:/usr/local/lib

export MERIC_NUM_THREADS=0
export MERIC_OUTPUT_DIR="ElmerData_meric"
export MERIC_DETAILED=1
export MERIC_AGGREGATE=1
export MERIC_MODE=2

export NPROCS=96
export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"

export SOLVER="cg"

# 2) Run
cd ${PROBLEM_PATH}
echo $PWD
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"

rm case.sif
cp ${ELMER_ROOT}/case_${SOLVER}.sif case.sif

stopHdeem
clearHdeem
startHdeem
sleep 1
stopHdeem

echo "Running Elmer (MERIC)"

clearHdeem
startHdeem

if [ ! -d ${PROBLEM_PATH}/winkel/partitioning.96 ]; then
    ${ELMER_ROOT}/gen_problem.sh
fi

echo "CF: ${MERIC_FREQUENCY}"
echo "UCF: ${MERIC_UNCORE_FREQUENCY}"
echo "SOLVER: ${SOLVER}"

# Core frequencies
for cf in $(seq 12 4 24) 25; do
    export MERIC_FREQUENCY="$cf"
    # Uncore frequencies
    for ucf in $(seq 12 4 28) 30; do
        export MERIC_UNCORE_FREQUENCY="$ucf"
        export MERIC_OUTPUT_FILENAME="${SOLVER}_${cf}_${ucf}"

        srun -n $NPROCS ${ELMER_HOME}/bin/ElmerSolver_mpi
    done
done

stopHdeem
sleep 1
checkHdeem

echo "Running Elmer (MERIC) done"
 
