#!/bin/sh

#SBATCH -t 0-02:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ondrej.vysocky@vsb.cz

#SBATCH --nodes=1
#SBATCH --tasks-per-node=28
#SBATCH --cpus-per-task=1

#SBATCH --partition=broadwell 
#SBATCH --reservation=p_readex_56
#SBATCH -A p_readex
#SBATCH --exclusive	
#SBATCH --mem-per-cpu=2200M
#SBATCH --comment="no_monitoring"
#SBATCH -J "kripke"

cd ../build
. ../readex_env/set_env_rdd.source
. ../environment.sh

export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM

echo "running kripke for readex-dyn-detect"
srun -n 28 ./kripke $KRIPKE_COMMAND
echo "running kripke done"

echo "running readex-dyn-detect"
echo "phase region = $2"
#readex-dyn-detect -t $1 -p $2 -c $3 -v $4 -w $5 scorep-*/profile.cubex
readex-dyn-detect -p "Loop" -t 0.01 scorep-*/profile.cubex
echo
echo "running readex-dyn-detect done" 

#cp -R scorep-* ../RESULTS/
#cp readex_config.xml ../RESULTS/
