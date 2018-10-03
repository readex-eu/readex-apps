#!/bin/sh

#SBATCH --time=00:40:00
#SBATCH --nodes=1  
#SBATCH --ntasks=28
#SBATCH --tasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --partition=broadwell
#SBATCH --comment="no_monitoring"
#SBATCH --mem-per-cpu=2200M
#SBATCH -J "FILTERfoam"
#SBATCH --mail-user=ondrej.vysocky@vsb.cz   # email address
#SBATCH --mail-type=ALL
#SBATCH -A p_readex
#SBATCH --reservation=p_readex_56

APP='srun -n 28 simpleFoam -parallel'

export FM_DIR=$( dirname "${BASH_SOURCE[0]}" )
if [ "$FM_DIR" == "." ]
then
        export FM_DIR=$(pwd)
fi
export FM_DIR=$FM_DIR/..

cd $FM_DIR
source readex_env/set_env_saf.source
source scripts_$READEX_MACHINE/environment.sh

cd $FM_DIR/OpenFOAM-v1612+/
source $FM_DIR/OpenFOAM-v1612+/etc/bashrc
cd ../../motorBike28/

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_FILTERING_FILE=scorep.filt

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt

iteration=0
result=1
while [ $result != 0 ]; do
  iteration=$(($iteration +1))
  echo "result = "$result
  # run the application
  $APP
  echo "ITERATION done."
  $FM_DIR/scripts_$READEX_MACHINE/do_scorep_autofilter_single.sh 0.1
  result=$?
  echo "scorep_autofilter_singe done ($result)."
  mv scorep-* $iteration"_scorep"
  cp scorep.filt $iteration"_scorep.filt"
done
echo "end."
