#!/bin/sh

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2500M  
#SBATCH -A p_readex
#SBATCH -J "READEX_kripke"

cd ../build
. ../readex_env/set_env_rdd.source
. ../environment.sh

#cp ../RESULTS/scorep.filt .

export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM
#export SCOREP_FILTERING_FILE=scorep.filt

echo "running kripke for readex-dyn-detect"
srun -n 24 ./kripke $KRIPKE_COMMAND
echo "running kripke done"

echo "running readex-dyn-detect"
echo "phase region = $2"
#readex-dyn-detect -t $1 -p $2 -c $3 -v $4 -w $5 scorep-*/profile.cubex
readex-dyn-detect -p "Loop" -t 0.01 scorep-*/profile.cubex
echo
echo "running readex-dyn-detect done" 

#cp -R scorep-* ../RESULTS/
#cp readex_config.xml ../RESULTS/
