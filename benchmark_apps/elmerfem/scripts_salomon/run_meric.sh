#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_meric
#PBS -l select=4:ncpus=24:mpiprocs=24:accelerator=false
#PBS -l walltime=01:00:00

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
    LOC=$PBS_O_WORKDIR
else
    LOC=$(pwd)
fi

cd $LOC


module load mkl
module load CMake/3.9.1

. readex_env/set_env_meric.source
. ../paths.source

cd ${ELMER_ROOT}


export LD_LIBRARY_PATH+=:/projects/p_readex/it4i/mericGCC/lib:/usr/local/lib

export MERIC_NUM_THREADS=0
export MERIC_OUTPUT_DIR="ElmerData_meric"
export MERIC_DETAILED=1
export MERIC_AGGREGATE=1
export MERIC_MODE=2

export NPROCS=96
export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"


# 2) Run
cd ${PROBLEM_PATH}
echo $PWD
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"

rm case.sif
cp ${ELMER_ROOT}/case_cg.sif case.sif

#stopHdeem
#clearHdeem
#startHdeem
#sleep 1
#stopHdeem

echo "Running Elmer (MERIC)"

#clearHdeem
#startHdeem

if [ ! -d ${PROBLEM_PATH}/winkel/partitioning.96 ]; then
    ${ELMER_ROOT}/gen_problem.sh
fi

# Core frequencies
for cf in $(seq 12 4 24) 25; do
    export MERIC_FREQUENCY="$cf"
    # Uncore frequencies
    for ucf in $(seq 12 4 28) 30; do
        export MERIC_UNCORE_FREQUENCY="$ucf"
        export MERIC_OUTPUT_FILENAME="${SOLVER}_${cf}_${ucf}"

        mpirun -n $NPROCS ${ELMER_HOME}/bin/ElmerSolver_mpi
    done
done

#stopHdeem
#sleep 1
#checkHdeem

echo "Running Elmer (MERIC) done"
 
