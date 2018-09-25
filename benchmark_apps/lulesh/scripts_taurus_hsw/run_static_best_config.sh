#!/bin/sh
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2500M
#SBATCH --comment="cpufreqchown"
#SBATCH -J "lulesh_sacct"
#SBATCH -A p_readex
#SBATCH --reservation=READEX
#SBATCH --output=static_average_best_config.out
#SBATCH --error=static_average_best_config.out

cd ..
REPEAT_COUNT=5
module purge
source ./readex_env/set_env_plain.source 

#module load scorep-hdeem/2016-12-20-hdeem-2.2.20ms
NUM_CPUS=24
#export SCOREP_METRIC_PLUGINS="hdeem_plugin"
#export SCOREP_METRIC_HDEEM_PLUGIN=*
#export SCOREP_METRIC_HDEEM_PLUGIN_VERBOSE=WARN
#export SCOREP_METRIC_HDEEM_PLUGIN_CONNECTION=INBAND
#export SCOREP_METRIC_HDEEM_PLUGIN_TIMER=BMC

change_frequency() {
for ((i = 0; i<$NUM_CPUS; i++))
do 
	cpufreq-set -c $i -f $1GHz
done
}

check_uncore_frequency() {
	x86a_read -n -i Intel_UNCORE_MIN_RATIO
   	x86a_read -n -i Intel_UNCORE_MAX_RATIO
    x86a_read -n -i Intel_UNCORE_CURRENT_RATIO
}




#export SCOREP_ENABLE_TRACING=true
#export SCOREP_ENABLE_PROFILING=false
#export SCOREP_TOTAL_MEMORY=3G
#export SCOREP_MPI_ENABLE_GROUPS=EXT

change_frequency 2.5 
cpufreq-set -c 16 -f 2.5GHz
cpufreq-set -c 17 -f 2.5GHz 
cpufreq-set -c 18 -f 2.5GHz 
cpufreq-set -c 19 -f 2.5GHz 
cpufreq-set -c 20 -f 2.5GHz 
cpufreq-set -c 21 -f 2.5GHz 
cpufreq-set -c 22 -f 2.5GHz 
cpufreq-set -c 23 -f 2.5GHz 
cpufreq-set -c 24 -f 2.5GHz 


x86a_write -n -c 0 -i Intel_UNCORE_MAX_RATIO -V 22
x86a_write -n -c 1 -i Intel_UNCORE_MAX_RATIO -V 22
x86a_write -n -c 0 -i Intel_UNCORE_MIN_RATIO -V 22 
x86a_write -n -c 1 -i Intel_UNCORE_MIN_RATIO -V 22

cpufreq-info
check_uncore_frequency

i=1
rm -rf TUNED_*
while [ $i -le $REPEAT_COUNT ]; do
  mkdir TUNED_$i
  export MEASURE_RAPL_TARGET="TUNED_$i"
  srun --ntasks 1 --ntasks-per-node 1 --cpus-per-task 24 measure-rapl ./lulesh2.0_plain -i 100 -s 150 
  i=$(echo "$i + 1" | bc)
done

i=1
total_time_rrl=0
total_energy_rrl=0
total_cpu_energy_rrl=0
while [ $i -lt $REPEAT_COUNT ]; do
   #echo "command sacct -j $SLURM_JOBID.$((i)) --format="JobID,CPUTimeRAW,ConsumedEnergyRaw""
  times_energys=$(sacct -j $SLURM_JOBID.$i --format="JobID,CPUTimeRAW,ConsumedEnergyRaw")
  i=$(echo "$i + 1" | bc)
  times_energys_array=(${times_energys[@]})
  time_step=${times_energys_array[7]}
  energy_step=${times_energys_array[8]}
  total_time_rrl=$(echo "${total_time_rrl} + ${time_step}" | bc)
  total_energy_rrl=$(echo "${total_energy_rrl} + ${energy_step}" | bc)
  echo "Time for job ${time_step}"
  echo "Energy for job ${energy_step}"
  for file in TUNED_$i/*
  do
    values=$( tail -1 $file | awk -F'[ ,]' '{print int($1)" "int($2)}' )
    values=(${values[@]})
    total_cpu_energy_rrl=$[ total_cpu_energy_rrl + ${values[0]} + ${values[1]} ]
  done
done

#echo "Total Plain Time = $total_time_plain, Total Plain Energy = $total_energy_plain"
#avg_time_plain=$(echo "$total_time_plain / $((REPEAT_COUNT-1))" | bc)
#avg_energy_plain=$(echo "$total_energy_plain / $((REPEAT_COUNT-1))" | bc)

echo "Total Static Time = $total_time_rrl, Total Static Energy = $total_energy_rrl"
avg_time_rrl=$(echo "$total_time_rrl/$((REPEAT_COUNT -1))" | bc)
avg_energy_rrl=$(echo "$total_energy_rrl/$((REPEAT_COUNT -1))" | bc)

echo "Average RRL Time = $avg_time_rrl"
echo "Average RRL Energy = $avg_energy_rrl"

rm -rf TUNED_*



