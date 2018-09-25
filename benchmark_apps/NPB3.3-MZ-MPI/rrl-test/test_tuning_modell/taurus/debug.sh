#!/bin/bash

module use /projects/p_readex/modules/

module load scorep/TRY_READEX_online_access_call_tree_extensions_r11444_bullxmpi_gcc5.3.0
module load readex-rrl/rrl-2016-11-08-gcc5.3.0-bullxmpi
module load pcp/pcp_2016-11-08-gcc3.0

echo $(hostname)

module list 

echo $(hostname)

CUR_DIR=`pwd`

export OMP_NUM_THREADS=24

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_TOTAL_MEMORY="1000M"

export SCOREP_METRIC_PLUGINS="hdeem_sync_plugin"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_VERBOSE="WARN"
export SCOREP_METRIC_HDEEM_SYNC_PLUGIN_STATS_TIMEOUT_MS="10"

export SCOREP_SUBSTRATE_PLUGINS='rrl'

export SCOREP_RRL_VERBOSE="DEBUG"
export SCOREP_TUNING_PLUGINS='OpenMPTP,cpu_freq_plugin'
export SCOREP_TUNING_CPU_FREQ_PLUGIN_VERBOSE="DEBUG"
export SCOREP_TUNING_OPENMPTP_PLUGIN_VERBOSE="DEBUG"
export SCOREP_RRL_TMM_PATH="$CUR_DIR/tuning_model.json"

cd ../../../bin/

gdb ./bt-mz.S.1

