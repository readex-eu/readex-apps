#!/bin/sh

#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2500M
#SBATCH -J "m_rrl_r"
#SBATCH -A p_readex
#SBATCH --output=bt-mz_sacct.out
#SBATCH --error=bt-mz_sacct.out

#####
# set to number of times to repeat experiments
REPEAT_COUNT=5
#####

module purge
source ./readex_env/set_env_plain.source

export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS=""
export SCOREP_RRL_PLUGINS=""
export SCOREP_RRL_TMM_PATH=""
export SCOREP_MPI_ENABLE_GROUPS=ENV

i=1
while [ $i -le $REPEAT_COUNT ]; do
  srun --cpu_bind=verbose,sockets ./bin/bt-mz.C.2_plain
  i=$(echo "$i + 1" | bc)
done

module purge
source ./readex_env/set_env_rrl.source

export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin,OpenMPTP"
export SCOREP_RRL_TMM_PATH="tuning_model.json"
export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_RRL_CHECK_IF_RESET="reset"
export SCOREP_EXPERIMENT_DIRECTORY=$(pwd)/scorep_m_rrl_r

i=1
while [ $i -le $REPEAT_COUNT ]; do
  srun --cpu_bind=verbose,sockets ./bin/bt-mz.C.2_ptf  
  i=$(echo "$i + 1" | bc)
done

i=1
total_time_plain=0
total_energy_plain=0
while [ $i -lt $REPEAT_COUNT ]; do
  times_energys=$(sacct -j $SLURM_JOBID.$i --format="JobID,CPUTimeRAW,ConsumedEnergyRaw")
  i=$(echo "$i + 1" | bc)
  times_energys_array=(${times_energys[@]})
  time_step=${times_energys_array[7]}
  energy_step=${times_energys_array[8]}
  total_time_plain=$(echo "${total_time_plain} + ${time_step}" | bc)
  total_energy_plain=$(echo "${total_energy_plain} + ${energy_step}" | bc)
done

i=1
total_time_rrl=0
total_energy_rrl=0
while [ $i -lt $REPEAT_COUNT ]; do
  times_energys=$(sacct -j $SLURM_JOBID.$((i+REPEAT_COUNT)) --format="JobID,CPUTimeRAW,ConsumedEnergyRaw")
  i=$(echo "$i + 1" | bc)
  times_energys_array=(${times_energys[@]})
  time_step=${times_energys_array[7]}
  energy_step=${times_energys_array[8]}
  total_time_rrl=$(echo "${total_time_rrl} + ${time_step}" | bc)
  total_energy_rrl=$(echo "${total_energy_rrl} + ${energy_step}" | bc)
done

avg_time_plain=$(echo "$total_time_plain / $((REPEAT_COUNT-1))" | bc)
avg_energy_plain=$(echo "$total_energy_plain / $((REPEAT_COUNT-1))" | bc)
avg_time_rrl=$(echo "$total_time_rrl / $((REPEAT_COUNT-1))" | bc)
avg_energy_rrl=$(echo "$total_energy_rrl / $((REPEAT_COUNT-1))" | bc)

echo ""
echo "Average Plain Time = $avg_time_plain, Average RRL Time = $avg_time_rrl"
echo "Average Plain Energy = $avg_energy_plain, Average RRL Energy = $avg_energy_rrl"

