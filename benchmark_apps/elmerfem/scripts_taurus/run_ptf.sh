#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --nodes=5
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH -p haswell
#SBATCH --mem-per-cpu=600M   # memory per CPU core
#SBATCH -J "READEX-run_ptf"
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH --comment="no_monitoring"
#SBATCH --mail-type=ALL
#SBATCH --reservation=READEX

module load mkl/2017

SCRIPT_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
. ${SCRIPT_DIR}/readex_env/set_env_ptf.source
. ${SCRIPT_DIR}/../paths.source

cd ${ELMER_ROOT}

#export SCOREP_SUBSTRATE_PLUGINS=PrintRegions

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export NPROCS=96
export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"

export SCOREP_SUBSTRATE_PLUGINS=rrl
export SCOREP_RRL_PLUGINS=cpu_freq_plugin,uncore_freq_plugin
export SCOREP_RRL_VERBOSE="WARN"

export SCOREP_METRIC_PLUGINS=hdeem_sync_plugin
export SCOREP_METRIC_PLUGINS_SEP=";"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN="*/E"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_CONNECTION="INBAND"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000
clearHdeem

#export SCOREP_METRIC_PLUGINS=x86_energy_sync_plugin
#export SCOREP_METRIC_PLUGINS_SEP=";"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN=*/E
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_CONNECTION="INBAND"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_VERBOSE="WARN"
#export SCOREP_METRIC_X86_ENERGY_SYNC_PLUGIN_STATS_TIMEOUT_MS=1000

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
