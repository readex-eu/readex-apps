#!/bin/sh

#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --partition=haswell
#SBATCH --reservation=READEX
#SBATCH -A p_readex
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2500M
#SBATCH --comment="no_monitoring"
#SBATCH -J "RDDfoam"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ondrej.vysocky@vsb.cz

APP='srun -n 24 simpleFoam -parallel'
PHASE_REG_NAME="iteration"

export FM_DIR=$( dirname "${BASH_SOURCE[0]}" )
if [ "$FM_DIR" == "." ]
then
        export FM_DIR=$(pwd)
fi
export FM_DIR=$FM_DIR/..

cd $FM_DIR
source readex_env/set_env_rdd.source
source scripts_$READEX_MACHINE/environment.sh

cd $FM_DIR/OpenFOAM-v1612+/
source $FM_DIR/OpenFOAM-v1612+/etc/bashrc
cd ../../motorBike24/


rm -rf scorep-*

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true
export SCOREP_MPI_ENABLE_GROUPS=ENV

echo "running simpleFoam for readex-dyn-detect"
# run the application
$APP
echo "running simpleFoam done"

echo "running readex-dyn-detect"
readex-dyn-detect -p $PHASE_REG_NAME -t 0.1  scorep-*/profile.cubex
echo
echo "running readex-dyn-detect done"
