#!/bin/sh

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=2  # number of processor cores (i.e. tasks)
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=6
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -J "e_rrl"   # job name
#SBATCH -A p_readex


module purge
source ./readex_env/set_env_ptf_hdeem.source
#module use /projects/p_readex/modules
#module load readex/ci_readex_bullxmpi1.2.8.4_gcc5.3.0

REPEAT_COUNT=5

INPUT_FILE=in3.data

export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS=""
export SCOREP_RRL_PLUGINS=""
export SCOREP_RRL_TMM_PATH=""
export SCOREP_MPI_ENABLE_GROUPS=ENV

i=1
while [ $i -le $REPEAT_COUNT ]; do
  srun ./miniMD_openmpi_plain -i $INPUT_FILE
  i=$(echo "$i + 1" | bc)
done

export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin,OpenMPTP"
export SCOREP_RRL_TMM_PATH="tuning_model.json"
export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_RRL_CHECK_IF_RESET="reset"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export SCOREP_EXPERIMENT_DIRECTORY=$(pwd)/scorep_e_rrl

i=1
while [ $i -le $REPEAT_COUNT ]; do
  srun ./miniMD_openmpi_ptf -i $INPUT_FILE
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
  #echo "${time_step}, ${energy_step}"
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
  #echo "${time_step}, ${energy_step}"
done

echo "Total Plain Time = $total_time_plain, Total Plain Energy = $total_energy_plain"
echo "Total RRL Time = $total_time_rrl, Total RRL Energy = $total_energy_rrl"

avg_time_plain=$(echo "$total_time_plain / $((REPEAT_COUNT-1))" | bc)
avg_energy_plain=$(echo "$total_energy_plain / $((REPEAT_COUNT-1))" | bc)
avg_time_rrl=$(echo "$total_time_rrl / $((REPEAT_COUNT-1))" | bc)
avg_energy_rrl=$(echo "$total_energy_rrl / $((REPEAT_COUNT-1))" | bc)

echo ""
echo "Average Plain Time = $avg_time_plain, Average RRL Time = $avg_time_rrl"
echo "Average Plain Energy = $avg_energy_plain, Average RRL Energy = $avg_energy_rrl"
