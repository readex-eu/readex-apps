#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_ptf
#PBS -l select=5:ncpus=24:mpiprocs=24:accelerator=false
#PBS -l walltime=02:00:00

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
    LOC=$PBS_O_WORKDIR
else
    LOC=$(pwd)
fi

cd $LOC

module load mkl
module load CMake/3.9.1

. readex_env/set_env_ptf.source
. ../paths.source

cd ${ELMER_ROOT}

#export SCOREP_SUBSTRATE_PLUGINS=PrintRegions

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export NPROCS=96
export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin
export SCOREP_RRL_VERBOSE="WARN"

export SCOREP_METRIC_PLUGINS=x86_energy_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN=*/E
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_FILTERING_FILE=scorep.filt

cp readex_config_extended.xml ${PROBLEM_PATH}
cd  ${PROBLEM_PATH}

cp ${ELMER_ROOT}/case_cg.sif case.sif
cp ${ELMER_ROOT}/scorep.filt .

echo "Generating Grid for Elmer..."
if [ ! -d ${PROBLEM_PATH}/winkel/partitioning.96 ]; then
    echo "Problem winkel generated."
    ${ELMER_ROOT}/gen_problem.sh
else
    echo "Problem winkel already exists."
fi

# run the application
echo "Running psc_frontend..."
psc_frontend --apprun="ElmerSolver_mpi" --mpinumprocs=${NPROCS} --ompnumthreads=1 --phase="ElmerSolver_core" --info=2 --config-file=readex_config_extended.xml --tune=readex_intraphase  --force-localhost --selective-info=AutotunePlugins,FrontendStateMachines 2>&1 | tee -a log_JOB.txt

echo "Running psc_frontend done."

cp rts.xml ${ELMER_ROOT}
cp tuning_model.json ${ELMER_ROOT}
