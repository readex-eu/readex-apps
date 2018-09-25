#!/bin/bash

module use /projects/p_readex/modules/

module load scorep/TRY_READEX_online_access_call_tree_extensions_r11444_bullxmpi_gcc5.3.0
module load readex-rrl/rrl-2016-11-08-gcc5.3.0-bullxmpi
module load pcp/pcp_2016-11-08-gcc3.0

echo $(hostname)

export OMP_NUM_THREADS=24

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_TOTAL_MEMORY="1000M"

export SCOREP_SUBSTRATE_PLUGINS='rrl'

export SCOREP_RRL_VERBOSE="TRACE"

export SCOREP_TUNING_PLUGINS='OpenMPTP,cpu_freq_plugin'
export SCOREP_TUNING_CPU_FREQ_PLUGIN_VERBOSE="DEBUG"
export SCOREP_TUNING_OPENMPTP_PLUGIN_VERBOSE="DEBUG"

export SCOREP_ONLINEACCESS_ENABLE=true

export SCOREP_ONLINEACCESS_REG_PORT=50100
export SCOREP_ONLINEACCESS_REG_HOST='localhost'
export SCOREP_ONLINEACCESS_BASE_PORT=50010
export SCOREP_ONLINEACCESS_APPL_NAME='appl'

scorep-online-access-registry 50100 test=./scenario_tuning_variable_fortran_foo_new_no_hdeem > online_access_output & scorep_oa=$!

CUR_DIR=`pwd`

cd ../../../bin/
srun -o "$CUR_DIR/tool_output" ./bt-mz.S.1
EXIT_CODE=$?


if [ -n "$(ps -p $scorep_oa -o pid=)" ]
then
	kill $scorep_oa
fi


exit $EXIT_CODE
