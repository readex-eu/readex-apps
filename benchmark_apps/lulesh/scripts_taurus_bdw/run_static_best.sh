#!/bin/sh
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --mem-per-cpu=2200M
#SBATCH --comment="cpufreqchown"
#SBATCH -J "lulesh_sacct"
#SBATCH -A p_readex
#SBATCH --reservation=p_readex_56
#SBATCH --comment="no_monitoring"
#SBATCH --output=static_average_1_nodes_new_1.out
#SBATCH --error=static_average_1_nodes_new_1.out

cd ..
REPEAT_COUNT=3
module purge
#source ./readex_env/set_env_plain.source 
source ./readex_env/set_env_saf.source

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/global/libraries/cpufrequtils/gcc5.3.0/lib/
#module load scorep-hdeem/2016-12-20-hdeem-2.2.20ms
NUM_CPUS=28
#export SCOREP_METRIC_PLUGINS="hdeem_plugin"
#export SCOREP_METRIC_HDEEM_PLUGIN=*
#export SCOREP_METRIC_HDEEM_PLUGIN_VERBOSE=WARN
#export SCOREP_METRIC_HDEEM_PLUGIN_CONNECTION=INBAND
#export SCOREP_METRIC_HDEEM_PLUGIN_TIMER=BMC

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


#change_frequency 2.4
#x86a_write -n -c 0 -i Intel_UNCORE_MAX_RATIO -V 27
#x86a_write -n -c 1 -i Intel_UNCORE_MAX_RATIO -V 27
#x86a_write -n -c 0 -i Intel_UNCORE_MIN_RATIO -V 27
#x86a_write -n -c 1 -i Intel_UNCORE_MIN_RATIO -V 27

i=1
rm -rf PLAIN_*
while [ $i -le $REPEAT_COUNT ]; do
  mkdir PLAIN_$i
  export MEASURE_RAPL_TARGET="PLAIN_$i"
  #srun --cpu_bind=verbose,sockets measure-rapl ./test/amg2013_plain -P 2 2 2 -r 40 40 40
  #srun --ntasks 8 --ntasks-per-node 1 --cpus-per-task 24 ./lulesh2.0_plain -i 500 -s 75
  srun measure-rapl ./lulesh2.0_plain -i 100 -s 150 
 i=$(echo "$i + 1" | bc)
done

#export SCOREP_ENABLE_TRACING=true
#export SCOREP_ENABLE_PROFILING=false
#export SCOREP_TOTAL_MEMORY=3G
#export SCOREP_MPI_ENABLE_GROUPS=EXT


#cpu_freq_list=(1.2 2.0 2.4 2.5)
#uncore_freq_list=(14 22 26 30)
#for i in "${cpu_freq_list[@]}"
#do
#	change_frequency $i
#	for j in "${uncore_freq_list[@]}"
#	do 
		#export MEASURE_RAPL_TARGET="TUNED_$sum"
#		cpufreq-info
#		x86a_write -n -c 0 -i Intel_UNCORE_MAX_RATIO -V $j
#		x86a_write -n -c 1 -i Intel_UNCORE_MAX_RATIO -V $j
#		x86a_write -n -c 0 -i Intel_UNCORE_MIN_RATIO -V $j
#		x86a_write -n -c 1 -i Intel_UNCORE_MIN_RATIO -V $j
#		check_uncore_frequency
		#srun measure-rapl ./lulesh2.0_plain -i 250 -s 75
		#srun -n 1 -c 24 --exclusive --mem-per-cpu 2500M -p haswell --reservation=READEX ./lulesh2.0_plain -i 50 -s 75 
		#sum=$sum + 1
#		mpiexec --np 1 --npernode 1 --cpus-per-proc 24 ./lulesh2.0_saf -i 350 -s 75
		#((sum++))
		#echo $sum
#	done
#done 

i=1
total_time_plain=0
total_energy_plain=0
total_cpu_energy_plain=0
while [ $i -lt $REPEAT_COUNT ]; do
#  echo "command sacct -j $SLURM_JOBID.$i --format="JobID,CPUTimeRAW,ConsumedEnergyRaw""
  times_energys=$(sacct -j $SLURM_JOBID.$i --format="JobID,CPUTimeRAW,ConsumedEnergyRaw")
  i=$(echo "$i + 1" | bc)
  times_energys_array=(${times_energys[@]})
  time_step=${times_energys_array[7]}
  energy_step=${times_energys_array[8]}
  echo "Job Time: $time_step"
  echo "Job Energy: $energy_step"
  total_time_plain=$(echo "${total_time_plain} + ${time_step}" | bc)
  total_energy_plain=$(echo "${total_energy_plain} + ${energy_step}" | bc)
  for file in PLAIN_$i/*
  do
    values=$( tail -1 $file | awk -F'[ ,]' '{print int($1)" "int($2)}' )
    values=(${values[@]})
    total_cpu_energy_plain=$[ total_cpu_energy_plain + ${values[0]} + ${values[1]} ]
  done
done

#i=1
#total_time_rrl=0
#total_energy_rrl=0
#total_cpu_energy_rrl=0
#while [ $i -lt $sum ]; do
#  echo "command sacct -j $SLURM_JOBID.$((i)) --format="JobID,CPUTimeRAW,ConsumedEnergyRaw""
#  times_energys=$(sacct -j $SLURM_JOBID.$((i)) --format="JobID,CPUTimeRAW,ConsumedEnergyRaw")
#  i=$(echo "$i + 1" | bc)
#  times_energys_array=(${times_energys[@]})
#  time_step=${times_energys_array[7]}
#  energy_step=${times_energys_array[8]}
  #total_time_rrl=$(echo "${total_time_rrl} + ${time_step}" | bc)
  #total_energy_rrl=$(echo "${total_energy_rrl} + ${energy_step}" | bc)
#  echo "Time for job ${time_step}"
#  echo "Energy for job ${energy_step}"
  #for file in TUNED_$i/*
  #do
   # values=$( tail -1 $file | awk -F'[ ,]' '{print int($1)" "int($2)}' )
    #values=(${values[@]})
    #total_cpu_energy_rrl=$[ total_cpu_energy_rrl + ${values[0]} + ${values[1]} ]
    #total_cpu_energy_rrl=${values[0]} + ${values[1]}
	#echo "Total Cpu Energy for one configuration : ${total_cpu_energy_rrl}"
  #done
#done


echo "Total Plain Time = $total_time_plain, Total Plain Energy = $total_energy_plain"
 
avg_time_plain=$(echo "$total_time_plain/$((REPEAT_COUNT-1))" | bc)
avg_energy_plain=$(echo "$total_energy_plain/$((REPEAT_COUNT-1))" | bc)

echo "Average Plain Time=$avg_time_plain"
echo "Average Plain Energy=$avg_energy_plain"

rm -rf PLAIN_*

#rm -rf PLAIN_*
#rm -rf TUNED_*


