#!/bin/bash

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "bem4i_saf"   # job name
#SBATCH --output=bem4i_saf.out
#SBATCH --error=bem4i_saf.out
###############################################################################

module purge
source ./readex_env/set_env_saf.source
source ./set_env_bem4i.source

export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

#SAF_t=0.001 #1 ms
#SAF_t=0.01 #10 ms
#SAF_t=0.1 #100 ms
SAF_t=1 #1000 ms

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  srun -N 1 -n 1 -c 24 ./dist/release_taurus/bem4i_saf 
  echo "Aplication run - done."
  ./do_scorep_autofilter_single.sh $SAF_t
  result=$?
  echo "scorep_autofilter_single done ($result)."
done

echo "end."
