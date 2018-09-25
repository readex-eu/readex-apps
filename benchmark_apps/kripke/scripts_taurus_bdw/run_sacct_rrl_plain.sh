#!/bin/sh

#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2200M
#SBATCH -J "kripke_sacct"
#SBATCH -A p_readex
#SBATCH --reservation=p_readex_56



#####
# set to number of times to repeat experiments
REPEAT_COUNT=5
#####

module purge
cd ../build
. ../readex_env/set_env_rrl.source
. ../environment.sh

export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS=""
export SCOREP_RRL_PLUGINS=""
export SCOREP_RRL_TMM_PATH=""
export SCOREP_MPI_ENABLE_GROUPS=ENV

i=1
rm -rf PLAIN_*
while [ $i -le $REPEAT_COUNT ]; do
  mkdir PLAIN_$i
  export MEASURE_RAPL_TARGET="PLAIN_$i"
  srun -n 28 measure-rapl ./kripke $KRIPKE_COMMAND
  i=$(echo "$i + 1" | bc)
done


export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin,OpenMPTP"
export SCOREP_RRL_TMM_PATH="tuning_model.json"
export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_RRL_CHECK_IF_RESET="reset"

i=1
rm -rf TUNED_*
while [ $i -le $REPEAT_COUNT ]; do
  mkdir TUNED_$i
  export MEASURE_RAPL_TARGET="TUNED_$i"
  srun -n 28 measure-rapl ./kripke $KRIPKE_COMMAND
  i=$(echo "$i + 1" | bc)
done

i=1
total_time_plain=0
total_energy_plain=0
total_cpu_energy_plain=0
while [ $i -lt $REPEAT_COUNT ]; do
  times_energys=$(sacct -j $SLURM_JOBID.$i --format="JobID,CPUTimeRAW,ConsumedEnergyRaw")
  i=$(echo "$i + 1" | bc)
  times_energys_array=(${times_energys[@]})
  time_step=${times_energys_array[7]}
  energy_step=${times_energys_array[8]}
  total_time_plain=$(echo "${total_time_plain} + ${time_step}" | bc)
  total_energy_plain=$(echo "${total_energy_plain} + ${energy_step}" | bc)
  for file in PLAIN_$i/*
  do
    values=$( tail -1 $file | awk -F'[ ,]' '{print int($1)" "int($2)}' )
    values=(${values[@]})
    total_cpu_energy_plain=$[ total_cpu_energy_plain + ${values[0]} + ${values[1]} ]
  done
done

i=1
total_time_rrl=0
total_energy_rrl=0
total_cpu_energy_rrl=0
while [ $i -lt $REPEAT_COUNT ]; do
  times_energys=$(sacct -j $SLURM_JOBID.$((i+REPEAT_COUNT)) --format="JobID,CPUTimeRAW,ConsumedEnergyRaw")
  i=$(echo "$i + 1" | bc)
  times_energys_array=(${times_energys[@]})
  time_step=${times_energys_array[7]}
  energy_step=${times_energys_array[8]}
  total_time_rrl=$(echo "${total_time_rrl} + ${time_step}" | bc)
  total_energy_rrl=$(echo "${total_energy_rrl} + ${energy_step}" | bc)
  for file in TUNED_$i/*
  do
    values=$( tail -1 $file | awk -F'[ ,]' '{print int($1)" "int($2)}' )
    values=(${values[@]})
    total_cpu_energy_rrl=$[ total_cpu_energy_rrl + ${values[0]} + ${values[1]} ]
  done
done

echo ""
echo "Total Plain Time = $total_time_plain, Total Plain Energy = $total_energy_plain"
echo "Total RRL Time = $total_time_rrl, Total RRL Energy = $total_energy_rrl"

avg_time_plain=$(echo "$total_time_plain / $((REPEAT_COUNT-1))" | bc)
avg_energy_plain=$(echo "$total_energy_plain / $((REPEAT_COUNT-1))" | bc)
avg_time_rrl=$(echo "$total_time_rrl / $((REPEAT_COUNT-1))" | bc)
avg_energy_rrl=$(echo "$total_energy_rrl / $((REPEAT_COUNT-1))" | bc)
avg_cpu_energy_plain=$(echo "$total_cpu_energy_plain / $((REPEAT_COUNT-1))" | bc)
avg_cpu_energy_rrl=$(echo "$total_cpu_energy_rrl / $((REPEAT_COUNT-1))" | bc)

echo ""
echo "Average Plain Time = $avg_time_plain, Average RRL Time = $avg_time_rrl"
echo "Average Plain Energy = $avg_energy_plain, Average RRL Energy = $avg_energy_rrl"
echo "Average Plain CPU Energy = $avg_cpu_energy_plain, Average RRL CPU Energy = $avg_cpu_energy_rrl"

rm -rf PLAIN_*
rm -rf TUNED_*

