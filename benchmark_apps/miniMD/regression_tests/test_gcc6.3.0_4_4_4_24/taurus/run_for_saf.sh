#!/bin/sh

source ../../init.sh
cd ${REL_PATH_APP_EXECUTION}

module purge
module use /projects/p_readex/modules
module load readex/ci_readex_bullxmpi1.2.8.4_gcc6.3.0

INPUT_FILE=in3.data
SAF_t=0.001

export SCOREP_FILTERING_FILE=scorep.filt
export SCOREP_TOTAL_MEMORY=3G

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  srun -N 1 -n 2 --ntasks-per-node=2 -c 12 --exclusive -p haswell --mem-per-cpu 2500M ./miniMD_openmpi_saf -i ${INPUT_FILE}
  echo "Aplication run - done."
  sh do_scorep_autofilter_single.sh $SAF_t
  result=$?
  echo "scorep_autofilter_singe done ($result)."
done

exit $result
