#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_saf
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

. readex_env/set_env_saf.source
. ../paths.source

cd ${ELMER_ROOT}



export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"

NPROCS=96
SAF_TIME=0.1

cd  ${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/
echo $PWD

if [ ! -d ${PROBLEM_PATH}/winkel/partitioning.96 ]; then
    echo "Problem winkel generated."
    ${ELMER_ROOT}/gen_problem.sh
else
    echo "Problem winkel already exists."
fi

rm case.sif
cp ${ELMER_ROOT}/case_cg.sif case.sif

export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm -f old_scorep*.filt
echo "" > scorep.filt

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  # run the application.. update this for different applications
  mpirun -n $NPROCS ${ELMER_HOME}/bin/ElmerSolver_mpi && echo "Elmer done"

  #cp -r scorep-* ${ELMER_ROOT}
  cp scorep.filt ${ELMER_ROOT}

  ${ELMER_ROOT}/do_scorep_autofilter_single.sh ${SAF_TIME}
  result=$?

  echo "scorep_autofilter_single done ($result)."
done
echo "end"

#cp scorep.filt ../RESULTS/
