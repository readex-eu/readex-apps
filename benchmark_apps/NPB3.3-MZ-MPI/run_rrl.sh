#!/bin/sh

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem-per-cpu=2500M
#SBATCH -J "bt-mz_rrl_r"
#SBATCH -A p_readex
#SBATCH --output=bt-mz_rrl.out
#SBATCH --error=bt-mz_rrl.out


energy_label="Energy"
srun -N 1 -n 1 --ntasks-per-node=1 -c 1 hostname >> bt-mz_rrl_host_names.out

module purge
source ./readex_env/set_env_plain.source

# start plain run
export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS=""
export SCOREP_RRL_PLUGINS=""
export SCOREP_RRL_TMM_PATH=""
export SCOREP_MPI_ENABLE_GROUPS=ENV

for run_id in $(seq 1 1); do
  srun -N 1 -n 1 --ntasks-per-node=1 -c 1 clearHdeem
  srun -N 1 -n 1 --ntasks-per-node=1 -c 1 startHdeem
  # run application
  start_time=$(($(date +%s%N)/1000000))
  srun --cpu_bind=verbose,sockets ./bin/bt-mz.C.2_plain
  stop_time=$(($(date +%s%N)/1000000))
  srun -N 1 -n 1 --ntasks-per-node=1 -c 1 stopHdeem
  srun -N 1 -n 1 --ntasks-per-node=1 -c 1 sleep 5
  exec < bt-mz_rrl_host_names.out
  while read host_name; do
    srun -N 1 -n 1 --ntasks-per-node=1 -c 1 --nodelist=$host_name checkHdeem >> bt-mz_rrl_temp_hdeem.out
  done
  
  energy_total=0
  if [ -e bt-mz_rrl_temp_hdeem.out ]; then
    exec < bt-mz_rrl_temp_hdeem.out
    while read max max_unit min min_unit average average_unit energy energy_unit; do
      if [ "$energy" == "$energy_label" ]; then
        read blade max_val min_val average_val energy_val
        energy_total=$(echo "${energy_total} + ${energy_val}" | bc)
      fi
    done
    time_total=$(echo "${stop_time} - ${start_time}" | bc)
    echo "bt-mz_rrl $run_id 1 2 2 12 $time_total $energy_total" >> ./bt-mz_rrl_plain_hdeem.out
    rm -rf bt-mz_rrl_temp_hdeem.out
  fi
  # end plain run
done



# start rrl-tuned run
module purge
source ./readex_env/set_env_rrl.source
export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin"
export SCOREP_RRL_TMM_PATH="./tuning_model.json"
export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_RRL_CHECK_IF_RESET="reset"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export SCOREP_TUNING_CPU_FREQ_PLUGIN_VERBOSE="DEBUG"
export SCOREP_TUNING_UNCORE_FREQ_PLUGIN_VERBOSE="DEBUG"
for run_id in $(seq 1 1); do
  date_name=$(date +%G%m%d%H%M%S%N)
  export SCOREP_EXPERIMENT_DIRECTORY="scorep-bt-mz_rrl_1_2_2_12_${date_name}"
  srun -N 1 -n 1 --ntasks-per-node=1 -c 1 clearHdeem
  srun -N 1 -n 1 --ntasks-per-node=1 -c 1 startHdeem
  # run application
  start_time=$(($(date +%s%N)/1000000))
  srun --cpu_bind=verbose,sockets ./bin/bt-mz.C.2_ptf
  stop_time=$(($(date +%s%N)/1000000))
  srun -N 1 -n 1 --ntasks-per-node=1 -c 1 stopHdeem
  srun -N 1 -n 1 --ntasks-per-node=1 -c 1 sleep 5
  exec < bt-mz_rrl_host_names.out
  while read host_name; do
    srun -N 1 -n 1 --ntasks-per-node=1 -c 1 --nodelist=$host_name checkHdeem >> bt-mz_rrl_temp_hdeem.out
  done
  
  energy_total=0
  if [ -e bt-mz_rrl_temp_hdeem.out ]; then
    exec < bt-mz_rrl_temp_hdeem.out
    while read max max_unit min min_unit average average_unit energy energy_unit; do
      if [ "$energy" == "$energy_label" ]; then
        read blade max_val min_val average_val energy_val
        energy_total=$(echo "${energy_total} + ${energy_val}" | bc)
      fi
    done
    time_total=$(echo "${stop_time} - ${start_time}" | bc)
    echo "bt-mz_rrl $run_id 1 2 2 12 $time_total $energy_total" >> ./bt-mz_rrl_rrl_hdeem.out
    rm -rf bt-mz_rrl_temp_hdeem.out
  fi
  # end rrl-tuned run
done

rm -rf bt-mz_rrl_host_names.out

