#!/bin/bash

#SBATCH -p haswell
#SBATCH --exclusive
#SBATCh --time=01:30:00
#SBATCH -c 24
#SBATCH -n 1
#SBATCh --cpu_bind=threads
#SBATCH --mem-per-cpu=2300

module use /projects/p_readex/modules/

module load scorep/scorep/ci_TRY_READEX_online_access_call_tree_extensions_bullxmpi1.2.8.4_gcc5.3.0
module load readex-rrl/ci_rrl_bullxmpi1.2.8.4_gcc5.3.0
module load pcp/ci_pcp_bullxmpi1.2.8.4_gcc5.3.0
module load scorep_plugins/control_plugins

echo $(hostname)

export OMP_NUM_THREADS=24

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_TOTAL_MEMORY="1000M"

export SCOREP_SUBSTRATE_PLUGINS='rrl'
export SCOREP_RRL_ATPS='PARAMETER1'
export SCOREP_METRIC_PLUGINS='scorep_substrate_rrl'
export SCOREP_METRIC_SCOREP_SUBSTRATE_RRL='*,ATP/PARAMETER1'
export SCOREP_RRL_VERBOSE="TRACE"

export SCOREP_TUNING_PLUGINS='OpenMPTP,cpu_freq_plugin'
export SCOREP_TUNING_CPU_FREQ_PLUGIN_VERBOSE="DEBUG"
export SCOREP_TUNING_OPENMPTP_PLUGIN_VERBOSE="DEBUG"

export SCOREP_ONLINEACCESS_ENABLE=true
export SCOREP_USER_ENABLE=true

export SCOREP_ONLINEACCESS_REG_PORT=50100
export SCOREP_ONLINEACCESS_REG_HOST='localhost'
export SCOREP_ONLINEACCESS_BASE_PORT=50010
export SCOREP_ONLINEACCESS_APPL_NAME='appl'

printenv |grep SCOREP_*

scorep-online-access-registry 50100 test=./scenario_tuning_function_c > online_access_output & scorep_oa=$!

CUR_DIR=`pwd`

gdb ./tuning_function

if [ -n "$(ps -p $scorep_oa -o pid=)" ]
then
	kill $scorep_oa
fi



