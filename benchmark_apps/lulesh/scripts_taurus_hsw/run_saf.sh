#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -J "lulesh_saf"   # job name
#SBATCH -A p_readex
#SBATCH --output=lulesh_saf.out
#SBATCH --error=lulesh_saf.out 
#SBATCH --reservation=READEX 

LOC=$(pwd)
module purge
source ./readex_env/set_env_saf.source

export SCOREP_FILTERING_FILE=scorep.filt
export OMP_NUM_THREADS=24

pathToFilter=$(pwd)

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

if [ "$READEX_INTEL" == "1" ]; then
	rm -rf old_scorep_icc.filt
	echo "" > scorep_icc.filt
fi

#SAF_t=0.001 #1 ms
#SAF_t=0.01 #10 ms
SAF_t=0.1 #100 ms

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  srun ./lulesh2.0_saf -p -i 100 -s 150
  echo "Aplication run - done."
  $LOC/do_scorep_autofilter_single.sh $SAF_t
  result=$?

  if [ "$READEX_INTEL" == "1" ] && [ $result != 0 ]; then
	  export FILTER_INTEL="-tcollect-filter="${pathToFilter}"/scorep_icc.filt"
	  ./compile_for_saf.sh 
  fi
  echo "scorep_autofilter_single done ($result)."
done

