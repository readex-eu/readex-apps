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
#SBATCH -J "bem4i_meric_opt"   # job name
#SBATCH --output=bem4i_meric_opt.out
#SBATCH --error=bem4i_meric_opt.out
###############################################################################

phase_name=main_phase
test_name=bem4i_meric
plain_test_out_file=${test_name}_plain_hdeem.out
meric_test_out_file=${test_name}_meric_hdeem.out
tuning_model_file='./energy.opts'
REPS=2
#NUMA=numactl" "--interleave=all"
NUMA=
PRG=srun" "-N" "1" "-n" "1" "-c" "24" "$NUMA" "./dist/release_taurus/bem4i_meric

# LOAD MODULES
module purge
source ./readex_env/set_env_meric.source
source ./set_env_bem4i.source

# SET MERIC OUTPUT
#export MERIC_COUNTERS=perfevent
export MERIC_MODE=0
export MERIC_DEBUG=0
export MERIC_CONTINUAL=1
export MERIC_AGGREGATE=1
export MERIC_DETAILED=0
export MERIC_NUM_THREADS=0
export MERIC_FREQUENCY=0
export MERIC_UNCORE_FREQUENCY=0

export MERIC_OUTPUT_DIR="meric_measurement"
rm -rf $MERIC_OUTPUT_DIR

# DEFAULT RUN
unset MERIC_REGION_OPTIONS
export MERIC_OUTPUT_FILENAME=default

for repeat in $(seq 1 $REPS)
do
  $PRG
done

# TUNED RUN
export MERIC_REGION_OPTIONS=$tuning_model_file
export MERIC_OUTPUT_FILENAME=tuned

for repeat in $(seq 1 $REPS)
do
  $PRG
done

echo 'default consumption [J]' > $plain_test_out_file
sed -n 's/Average energy consumption \[J\],\(.*\)/\1/p' ${MERIC_OUTPUT_DIR}/${phase_name}/*_default.csv >> $plain_test_out_file 
echo 'default runtime [s]' >> $plain_test_out_file
sed -n 's/Maximal runtime of function \[s\],\(.*\)/\1/p' ${MERIC_OUTPUT_DIR}/${phase_name}/*_default.csv >> $plain_test_out_file

echo 'tuned consumption [J]' > $meric_test_out_file
sed -n 's/Average energy consumption \[J\],\(.*\)/\1/p' ${MERIC_OUTPUT_DIR}/${phase_name}/*_tuned.csv >> $meric_test_out_file 
echo 'tuned runtime [s]' >> $meric_test_out_file
sed -n 's/Maximal runtime of function \[s\],\(.*\)/\1/p' ${MERIC_OUTPUT_DIR}/${phase_name}/*_tuned.csv >> $meric_test_out_file

rm -rf $MERIC_OUTPUT_DIR
rm -rf ${MERIC_OUTPUT_DIR}Counters
