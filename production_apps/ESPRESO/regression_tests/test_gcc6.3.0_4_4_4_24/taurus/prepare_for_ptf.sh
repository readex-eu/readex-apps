#!/bin/sh

rm slurm*.out

source ../../init.sh
rm ${REL_TUNING_MODEL_FILE_NAME}
cd ${REL_PATH_APP_EXECUTION}

#module purge
#module use /projects/p_readex/modules
#module load readex/ci_readex_bullxmpi1.2.8.4_gcc6.3.0

###########
rdd_config_file="readex_config.xml"
ptf_config_file="ptf_readex_config.xml"
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

  echo "    </readex-dyn-detect>" >> $ptf_config_file
  echo "    <tuningParameter>" >> $ptf_config_file
  echo "        <frequency>" >> $ptf_config_file
  echo "            <min_freq>1200000</min_freq>" >> $ptf_config_file
  echo "            <max_freq>2400000</max_freq>" >> $ptf_config_file
  echo "            <freq_step>1200000</freq_step>" >> $ptf_config_file
  echo "        </frequency>" >> $ptf_config_file
  echo "        <uncore>" >> $ptf_config_file
  echo "            <min_freq>12</min_freq>" >> $ptf_config_file
  echo "            <max_freq>24</max_freq>" >> $ptf_config_file
  echo "            <freq_step>12</freq_step>" >> $ptf_config_file
  echo "        </uncore>" >> $ptf_config_file
  echo "        <openMPThreads>" >> $ptf_config_file
  echo "            <lower_value>24</lower_value>" >> $ptf_config_file
  echo "            <step>24</step>" >> $ptf_config_file
  echo "        </openMPThreads>" >> $ptf_config_file
  echo "    </tuningParameter>" >> $ptf_config_file
  echo "    <objectives>" >> $ptf_config_file
  echo "        <objective>Energy</objective>" >> $ptf_config_file
  echo "        <objective>Time</objective>" >> $ptf_config_file
  echo "        <objective>EDP</objective>" >> $ptf_config_file
  echo "        <objective>ED2P</objective>" >> $ptf_config_file
  echo "        <objective>CPUEnergy</objective>" >> $ptf_config_file
  echo "    </objectives>" >> $ptf_config_file
  echo "    <periscope>" >> $ptf_config_file
  echo "        <metrics>" >> $ptf_config_file
  echo "            <node_energy>hdeem/BLADE/E</node_energy>" >> $ptf_config_file
  echo "            <cpu0_energy>hdeem/CPU0/E</cpu0_energy>" >> $ptf_config_file
  echo "            <cpu1_energy>hdeem/CPU1/E</cpu1_energy>" >> $ptf_config_file
  echo "        </metrics>" >> $ptf_config_file
  echo "        <searchAlgorithm>" >> $ptf_config_file
  echo "            <name>exhaustive</name>" >> $ptf_config_file
  echo "        </searchAlgorithm>" >> $ptf_config_file
  echo "    </periscope>" >> $ptf_config_file
  echo "    <scorep>" >> $ptf_config_file
  echo "        <tuningSubstrate>rrl</tuningSubstrate>" >> $ptf_config_file
  echo "    </scorep>" >> $ptf_config_file
  echo "</Configuration>" >> $ptf_config_file
fi
###########
rm espreso_ptf
sh compile_for_ptf.sh
