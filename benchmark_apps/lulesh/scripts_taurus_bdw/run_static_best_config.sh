#!/bin/sh

#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --mem-per-cpu=2200M
#SBATCH --comment="cpufreqchown"
#SBATCH -J "lulesh_sacct"
#SBATCH -A p_readex
#SBATCH --comment="no_monitoring"
#SBATCH --reservation=p_readex_56
#SBATCH --output=static_one_node_plain_best_new.out
#SBATCH --error=static_one_node_plain_best_new.out

cd ..
REPEAT_COUNT=3
module purge
source ./readex_env/set_env_plain.source 


NUM_CPUS=28

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/global/libraries/cpufrequtils/gcc5.3.0/lib/ 

change_frequency() {
for ((i = 0; i<$NUM_CPUS; i++))
do 
	/sw/global/libraries/cpufrequtils/gcc5.3.0/bin/cpufreq-set -c $i -f $1GHz
done
}

check_uncore_frequency() {
	x86a_read -n -i Intel_UNCORE_MIN_RATIO
   	x86a_read -n -i Intel_UNCORE_MAX_RATIO
    x86a_read -n -i Intel_UNCORE_CURRENT_RATIO
}

change_frequency 2.4
x86a_write -n -c 0 -i Intel_UNCORE_MAX_RATIO -V 19
x86a_write -n -c 1 -i Intel_UNCORE_MAX_RATIO -V 19
x86a_write -n -c 0 -i Intel_UNCORE_MIN_RATIO -V 19
x86a_write -n -c 1 -i Intel_UNCORE_MIN_RATIO -V 19


/sw/global/libraries/cpufrequtils/gcc5.3.0/bin/cpufreq-info


i=1
#rm -rf TUNED_*
while [ $i -le $REPEAT_COUNT ]; do
  mkdir TUNED_$i
  #export MEASURE_RAPL_TARGET="TUNED_$i"
  #srun --cpu_bind=verbose,sockets measure-rapl ./lulesh2.0_ptf -i 200 -s 75
  srun -n 1 -c 28 measure-rapl ./lulesh2.0_plain -i 100 -s 150 
  i=$(echo "$i + 1" | bc)
done

i=1
total_time_rrl=0
total_energy_rrl=0
total_cpu_energy_rrl=0
while [ $i -lt $REPEAT_COUNT ]; do
  #echo "command: sacct -j $SLURM_JOBID.$((i)) --format="JobID,CPUTimeRAW,ConsumedEnergyRaw""
  times_energys=$(sacct -j $SLURM_JOBID.$i --format="JobID,CPUTimeRAW,ConsumedEnergyRaw")
  i=$(echo "$i + 1" | bc)
  times_energys_array=(${times_energys[@]})
  time_step=${times_energys_array[7]}
  energy_step=${times_energys_array[8]}
  total_time_rrl=$(echo "${total_time_rrl} + ${time_step}" | bc)
  total_energy_rrl=$(echo "${total_energy_rrl} + ${energy_step}" | bc)
  echo "Job time: ${time_step}"
  echo "Job Energy: ${energy_step}"
  #for file in TUNED_$i/*
  #do
  #  values=$( tail -1 $file | awk -F'[ ,]' '{print int($1)" "int($2)}' )
  #  values=(${values[@]})
  #  total_cpu_energy_rrl=$[ total_cpu_energy_rrl + ${values[0]} + ${values[1]} ]
  #done
done

#echo "Total Plain Time = $total_time_plain, Total Plain Energy = $total_energy_plain"
avg_time_rrl=$(echo "$total_time_rrl/$((REPEAT_COUNT-1))" | bc)
avg_energy_rrl=$(echo "$total_energy_rrl/$((REPEAT_COUNT-1))" | bc)

echo "Average RRL Time = $avg_time_rrl"
echo "Average RRL Energy = $avg_energy_rrl" 

