#!/bin/sh

# value "" if not used
CORE_MIN_FREQ="1200000"
CORE_MAX_FREQ="3000000"
CORE_FREQ_STEP="100000"

# value "" if not used
UNCORE_MIN_FREQ="10"
UNCORE_MAX_FREQ="30"
UNCORE_FREQ_STEP="1"

# value "" if not used
OMPT_LOWER_VALUE="" #"24"
OMPT_STEP="" #"24"

# if "" then default is "Energy"; remove objectives as required
OBJECTIVES="Energy Time EDP ED2P CPUEnergy"

# for HDEEM
METRIC_PLUGIN="hdeem_sync_plugin"
NODE_ENERGY="hdeem/BLADE/E"
CPU0_ENERGY="hdeem/CPU0/E"
CPU1_ENERGY="hdeem/CPU1/E"

# for RAPL
#METRIC_PLUGIN="x86_energy_sync_plugin"
#NODE_ENERGY="x86_energy/BLADE/E"
#CPU0_ENERGY="x86_energy/CORE0/E"
#CPU1_ENERGY="x86_energy/CORE1/E"

# if "" then default is "exhaustive"
SEARCH_ALGORITHM="exhaustive"
# SEARCH_ALGORITHM="random 2"
# SEARCH_ALGORITHM="individual 2"
# SEARCH_ALGORITHM="gde3 10 10 20"

SCOREP_TUNING_SUBSTRATE="rrl"

TUNING_MODEL_PATH="tuning_model.json"

###########
rdd_config_file="readex_config.xml" # readex_config file from readex-dyn-detect
ptf_config_file="ptf_readex_config.xml" # extended readex_config file for PTF
rm $ptf_config_file


if [ -e $rdd_config_file ]; then
  exec < $rdd_config_file
  IFS=''
  while read line; do
    if [ "$line" = "    </readex-dyn-detect>" ]; then
      break
    else
      echo $line >> $ptf_config_file
    fi
  done

IFS=' ' read -r -a objs <<< "$OBJECTIVES"
IFS=' ' read -r -a algo <<< "$SEARCH_ALGORITHM"

  echo "    </readex-dyn-detect>" >> $ptf_config_file
  echo "    <tuningParameter>" >> $ptf_config_file
  if [ "$CORE_MIN_FREQ" != "" ] && [ "$CORE_MAX_FREQ" != "" ] && [ "$CORE_FREQ_STEP" != "" ]; then
    echo "        <frequency>" >> $ptf_config_file
    echo "            <min_freq>$CORE_MIN_FREQ</min_freq>" >> $ptf_config_file
    echo "            <max_freq>$CORE_MAX_FREQ</max_freq>" >> $ptf_config_file
    echo "            <freq_step>$CORE_FREQ_STEP</freq_step>" >> $ptf_config_file
    echo "        </frequency>" >> $ptf_config_file
  fi
  if [ "$UNCORE_MIN_FREQ" != "" ] && [ "$UNCORE_MAX_FREQ" != "" ] && [ "$UNCORE_FREQ_STEP" != "" ]; then
    echo "        <uncore>" >> $ptf_config_file
    echo "            <min_freq>$UNCORE_MIN_FREQ</min_freq>" >> $ptf_config_file
    echo "            <max_freq>$UNCORE_MAX_FREQ</max_freq>" >> $ptf_config_file
    echo "            <freq_step>$UNCORE_FREQ_STEP</freq_step>" >> $ptf_config_file
    echo "        </uncore>" >> $ptf_config_file
  fi
  if [ "$OMPT_LOWER_VALUE" != "" ] && [ "$OMPT_STEP" != "" ]; then
    echo "        <openMPThreads>" >> $ptf_config_file
    echo "            <lower_value>$OMPT_LOWER_VALUE</lower_value>" >> $ptf_config_file
    echo "            <step>$OMPT_STEP</step>" >> $ptf_config_file
    echo "        </openMPThreads>" >> $ptf_config_file
  fi
  echo "    </tuningParameter>" >> $ptf_config_file
  echo "    <objectives>" >> $ptf_config_file
  if [ "$OBJECTIVES" = "" ]; then
    echo "        <objective>Energy</objective>" >> $ptf_config_file
  else
    for obj in "${objs[@]}"; do
      echo "        <objective>$obj</objective>" >> $ptf_config_file
    done
  fi
  echo "    </objectives>" >> $ptf_config_file
  echo "    <periscope>" >> $ptf_config_file
  echo "        <metricPlugin>" >> $ptf_config_file
  echo "            <name>$METRIC_PLUGIN</name>" >> $ptf_config_file
  echo "        </metricPlugin>" >> $ptf_config_file
  echo "        <metrics>" >> $ptf_config_file
  echo "            <node_energy>$NODE_ENERGY</node_energy>" >> $ptf_config_file
  echo "            <cpu0_energy>$CPU0_ENERGY</cpu0_energy>" >> $ptf_config_file
  echo "            <cpu1_energy>$CPU1_ENERGY</cpu1_energy>" >> $ptf_config_file
  echo "        </metrics>" >> $ptf_config_file
  echo "        <searchAlgorithm>" >> $ptf_config_file
  if [ "$SEARCH_ALGORITHM" = "" ]; then
    echo "            <name>exhaustive</name>" >> $ptf_config_file
  else
    if [ "${algo[0]}" = "exhaustive" ]; then
      echo "            <name>exhaustive</name>" >> $ptf_config_file
    elif [ "${algo[0]}" = "random" ]; then
      echo "            <name>random</name>" >> $ptf_config_file
      echo "                <samples>${algo[1]}</samples>" >> $ptf_config_file
    elif [ "${algo[0]}" = "individual" ]; then
      echo "            <name>individual</name>" >> $ptf_config_file
      echo "                <keep>${algo[1]}</keep>" >> $ptf_config_file
    elif [ "${algo[0]}" = "gde3" ]; then
      echo "            <name>gde3</name>" >> $ptf_config_file
      echo "                <popupationSize>${algo[1]}</populationSize>" >> $ptf_config_file
      echo "                <maxGeneration>${algo[2]}</maxGeneration>" >> $ptf_config_file
      echo "                <timer>${algo[3]}</timer>" >> $ptf_config_file
    else
      echo "            <name>exhaustive</name>" >> $ptf_config_file
    fi
  fi
  echo "        </searchAlgorithm>" >> $ptf_config_file
  echo "        <tuningModel>$TUNING_MODEL_PATH</tuningModel>" >> $ptf_config_file
  echo "    </periscope>" >> $ptf_config_file
  echo "    <scorep>" >> $ptf_config_file
  echo "        <tuningSubstrate>$SCOREP_TUNING_SUBSTRATE</tuningSubstrate>" >> $ptf_config_file
  echo "    </scorep>" >> $ptf_config_file
  echo "</Configuration>" >> $ptf_config_file
fi
###########
