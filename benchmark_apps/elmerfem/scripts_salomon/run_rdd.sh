#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_rdd
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
#module load PAPI

. readex_env/set_env_rdd.source
. ../paths.source

export SCOREP_PROFILING_FORMAT=cube_tuple
#export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
#export SCOREP_FILTERING_FILE=scorep.filt

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
mpirun -n ${NPROCS} ${ELMER_HOME}/bin/ElmerSolver_mpi
echo "Running  ElmerSolver done."

echo "Running readex-dyn-detect..."
readex-dyn-detect -p "ElmerSolver_core" -t "${RDD_TIME}" -c 10 -v 10 -w 10 scorep-*/profile.cubex
#readex-dyn-detect -t "${RDD_TIME}" -c 10 -v 10 -w 10 scorep-*/profile.cubex

echo "RDD result = $?"

echo "Running readex-dyn-detect done." 

#cp -R scorep-* ../RESULTS/
cp readex_config.xml ${ELMER_ROOT}
