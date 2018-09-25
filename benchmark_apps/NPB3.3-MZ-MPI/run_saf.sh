#!/bin/sh

##SBATCH --time=1:00:00   # walltime
##SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
##SBATCH --ntasks-per-node=2
##SBATCH --cpus-per-task=12
##SBATCH --exclusive
##SBATCH --partition=haswell
##SBATCH --comment="cpufreqchown"
##SBATCH --mem-per-cpu=2500M   # memory per CPU core
##SBATCH -A p_readex
##SBATCH -J ""$app"_saf"   # job name
##SBATCH --output="$app"_saf.out
##SBATCH --error="$app"_saf.out
###############################################################################


#source env/set_env_npb_ptf.sh
source readex_env/set_env_saf.source

export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

#SAF_t=0.001 #1 ms
SAF_t=0.01 #10 ms
#SAF_t=0.1 #100 ms
app=bt-mz.C.2_saf
result=1
while [ $result != 0 ]; do
  echo "result = "$result
srun -N 1 -n 2 -c 12 --exclusive -p interactive --mem-per-cpu 2500M bin/$app
  echo "Aplication run - done."
  ./do_scorep_autofilter_single.sh $SAF_t
  result=$?
  echo "scorep_autofilter_singe done ($result)."
done

echo "end."

