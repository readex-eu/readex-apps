#!/bin/sh

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
#SBATCH -J "e_saf"   # job name
#SBATCH --output=miniMD_saf.out
#SBATCH --error=miniMD_saf.err
###############################################################################

module purge
source ./readex_env/set_env_saf.source
INPUT_FILE=in3.data

export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

SAF_t=0.001 #1 ms
#SAF_t=0.01 #10 ms
#SAF_t=0.1 #100 ms

result=1
while [ $result != 0 ]; do
  echo "result = "$result
  srun ./miniMD_openmpi_saf -i $INPUT_FILE
  echo "Aplication run - done."
  sh do_scorep_autofilter_single.sh $SAF_t
  result=$?
  echo "scorep_autofilter_single done ($result)."
done

