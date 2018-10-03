#!/bin/bash 

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N espreso_meric
#PBS -l select=1:ncpus=24:mpiprocs=2:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=24:00:00  
#PBS -m be



hostname
if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
cd $FM_DIR/../../

source env/readex_env/set_env_meric.source
source env/modules.readex_salomon

source env/paths.default
source env/threading.default 12
################################################################################

# SET MERIC OUTPUT FORMAT
export MERIC_MODE=1
export MERIC_COUNTERS=perfevent
export MERIC_CONTINUAL=1
export MERIC_DETAILED=0
export MERIC_AGGREGATE=1
export MERIC_OUTPUT_DIR=$SCRATCH/DELETE

export MERIC_FREQUENCY=25
export MERIC_UNCORE_FREQUENCY=25
export MERIC_NUM_THREADS=12
export MERIC_OUTPUT_FILENAME=$MERIC_FREQUENCY"_"$MERIC_UNCORE_FREQUENCY"_CONFIG"


# RUN ONE TEST
mpirun -n 2 ./espreso_meric -c espreso.ecf

export MERIC_OUTPUT_DIR=$SCRATCH/ESPRESOgccFIN

# FOR EACH SETTINGS
for proc in 2
do
	for thread in {12..2..2}
	do
		for cpu_freq in 25 {24..12..2} 
		do
			for uncore_freq in {30..12..2} 
			do
				# OUTPUT NAMES
				export MERIC_OUTPUT_FILENAME=$thread"_"$cpu_freq"_"$uncore_freq"_CONFIG"

				# TEST SETTINGS
				export MERIC_FREQUENCY=$cpu_freq
				export MERIC_UNCORE_FREQUENCY=$uncore_freq
				export MERIC_NUM_THREADS=$thread

				echo "Output file: "  $MERIC_OUTPUT_FILENAME >> LOGmeric
				echo 
				
				mpirun -n 2 ./espreso_meric -c espreso.ecf | tee -a LOGmeric

			done
		done
	done
done



