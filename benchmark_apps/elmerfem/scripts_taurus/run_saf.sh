#!/bin/sh

#SBATCH --time=00:30:00   # walltime
#SBATCH --nodes=4  
#SBATCH --ntasks=24
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem=24000
#SBATCH -J "READEX-run_saf"   # job name
#SBATCH --mail-type=ALL
#SBATCH -A p_readex
#SBATCH --reservation=READEX

module load mkl/2017

SCRIPT_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
. ${SCRIPT_DIR}/readex_env/set_env_saf.source
. ${SCRIPT_DIR}/../paths.source

#cd ${ELMER_ROOT}

export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"

NPROCS=96
SAF_TIME=0.1

cd  ${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/

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
  srun -n $NPROCS ${ELMER_HOME}/bin/ElmerSolver_mpi && echo "Elmer done"

  #cp -r scorep-* ${ELMER_ROOT}
  cp scorep.filt ${ELMER_ROOT}

  ${ELMER_ROOT}/do_scorep_autofilter_single.sh ${SAF_TIME}
  result=$?

  echo "scorep_autofilter_single done ($result)."
done
echo "end"

#cp scorep.filt ../RESULTS/
