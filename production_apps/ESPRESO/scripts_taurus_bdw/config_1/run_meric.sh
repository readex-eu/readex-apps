#!/bin/bash 
#SBATCH -t 1-12:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ondrej.vysocky@vsb.cz

#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=14

#SBATCH --partition=broadwell 
#SBATCH --reservation=p_readex_44
#SBATCH -A p_readex	#to account your compute time on the readex project
#SBATCH --exclusive	
#SBATCH --mem-per-cpu=2200M
#SBATCH --comment="no_monitoring"

hostname
cd ../../

source env/readex_env/set_env_meric.source
source env/environment.sh

source env/paths.default
source env/threading.default 14
################################################################################

# SET MERIC OUTPUT FORMAT
export MERIC_MODE=1
#export MERIC_COUNTERS=perfevent
export MERIC_CONTINUAL=1
export MERIC_DETAILED=0
export MERIC_AGGREGATE=1
export MERIC_OUTPUT_DIR=$SCRATCH/DELETE

export MERIC_FREQUENCY=25
export MERIC_UNCORE_FREQUENCY=25
export MERIC_NUM_THREADS=12
export MERIC_OUTPUT_FILENAME=$MERIC_FREQUENCY"_"$MERIC_UNCORE_FREQUENCY"_CONFIG"


# RUN ONE TEST
srun -n 2 -c 14 ./espreso_meric -c espreso.ecf

export MERIC_OUTPUT_DIR=$SCRATCH/ESPRESOgccBDW

# FOR EACH SETTINGS
for proc in 2
do
	for thread in {14..2..4}
	do
		for cpu_freq in {24..12..2} 
		do
			for uncore_freq in 27 {26..12..2} 
			do
				# OUTPUT NAMES
				export MERIC_OUTPUT_FILENAME=$thread"_"$cpu_freq"_"$uncore_freq"_CONFIG"

				# TEST SETTINGS
				export MERIC_FREQUENCY=$cpu_freq
				export MERIC_UNCORE_FREQUENCY=$uncore_freq
				export MERIC_NUM_THREADS=$thread

				echo "Output file: "  $MERIC_OUTPUT_FILENAME
				echo 
				
				ret=1
				while [ "$ret" -ne 0 ]
				do
					srun -n 2 -c 14 ./espreso_meric -c espreso.ecf | tee -a LOGmeric
					ret=$?
				done

			done
		done
	done
done



