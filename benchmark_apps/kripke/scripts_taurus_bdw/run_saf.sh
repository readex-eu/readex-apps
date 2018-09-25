#!/bin/sh

#SBATCH -t 0-02:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ondrej.vysocky@vsb.cz

#SBATCH --nodes=1
#SBATCH --tasks-per-node=28
#SBATCH --cpus-per-task=1

#SBATCH --partition=broadwell 
#SBATCH --reservation=p_readex_56
#SBATCH -A p_readex
#SBATCH --exclusive	
#SBATCH --mem-per-cpu=2200M
#SBATCH --comment="no_monitoring"
#SBATCH -J "kripke"

cd ../build
. ../readex_env/set_env_saf.source
. ../environment.sh

export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

if [ "$READEX_INTEL" == "1" ]; then
  rm -rf old_scorep_icc.filt
  echo "" > scorep_icc.filt
fi


iter=0
result=1
while [ $result != 0 ]; do
  iter=$[ $iter +1 ]
  echo "result = "$result
  # run the application.. update this for different applications
  srun -n 28 ./kripke $KRIPKE_COMMAND
  echo "kripke done."
  ./do_scorep_autofilter_single.sh 0.001 #$1
  result=$?
#  cp scorep_icc.filt scorep_icc_"$iter".filt
#  cp scorep_gcc.filt scorep_gcc_"$iter".filt
  echo "scorep_autofilter_singe done ($result)."
done
echo "end." 

#cp scorep.filt ../RESULTS/
