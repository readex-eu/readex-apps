#!/bin/sh

source ../../init.sh
cd ${REL_PATH_APP_EXECUTION}

#module purge
#module use /projects/p_readex/modules
#module load readex/ci_readex_bullxmpi1.2.8.4_gcc6.3.0

source env/modules.taurus.atp
source env/paths.default
source env/threading.default 12

INPUT_FILE=espreso.ecf
SAF_t=0.01

export SCOREP_FILTERING_FILE=scorep.filt
export SCOREP_TOTAL_MEMORY=3G

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  srun -N 1 -n 2 --ntasks-per-node=2 -c 12 --exclusive -p haswell --mem-per-cpu 2500M ./espreso_saf
  echo "Aplication run - done."
  sh do_scorep_autofilter_single.sh $SAF_t
  result=$?
  echo "scorep_autofilter_singe done ($result)."
done

exit $result
