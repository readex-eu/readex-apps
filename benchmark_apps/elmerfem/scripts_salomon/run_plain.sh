#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_plain
#PBS -l select=4:ncpus=24:mpiprocs=24:accelerator=false
#PBS -l walltime=10:00:00

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
    LOC=$PBS_O_WORKDIR
else
    LOC=$(pwd)
fi

cd $LOC


module load mkl
module load CMake

. readex_env/set_env_plain.source
. ../paths.source

cd ${ELMER_ROOT}

export NPROCS=96

export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"

# 2) Run
cd ${PROBLEM_PATH}
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"
rm -f case.sif
cp ${ELMER_ROOT}/case_cg.sif case.sif


#stopHdeem
#clearHdeem
#startHdeem
#sleep 1
#stopHdeem


echo "Running Elmer (plain)"

#clearHdeem
#startHdeem

if [ ! -d ${PROBLEM_PATH}/winkel/partitioning.96 ]; then
    ${ELMER_ROOT}/gen_problem.sh
fi

mpirun -np $NPROCS ElmerSolver_mpi

#stopHdeem
#sleep 1
#checkHdeem

echo "Running Elmer (plain) done"

