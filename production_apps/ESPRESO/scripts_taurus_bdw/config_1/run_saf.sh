#!/bin/sh

#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ondrej.vysocky@vsb.cz

#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=14

#SBATCH --partition=broadwell 
#SBATCH --reservation=p_readex_56
#SBATCH -A p_readex	#to account your compute time on the readex project
#SBATCH --exclusive	
#SBATCH --mem-per-cpu=2200M
#SBATCH --comment="no_monitoring"

cd ../../
source env/readex_env/set_env_saf.source
source env/environment.sh

source env/paths.default
source env/threading.default 14


export SCOREP_TOTAL_MEMORY=3G
export SCOREP_FILTERING_FILE=scorep.filt

#SAF_t=0.001 #1 ms
#SAF_t=0.01 #10 ms
SAF_t=0.1 #100 ms

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

if [ "$READEX_INTEL" == "1" ]; then
  rm -rf old_scorep_icc.filt
  echo "" > scorep_icc.filt
fi

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  srun -n 2 -c 14 ./espreso_saf -c espreso_scaling.ecf
  echo "Aplication run - done."
  ./do_scorep_autofilter_single.sh $SAF_t
  result=$?
  echo "scorep_autofilter_singe done ($result)."
done

echo "end of scorep-autofilter."
