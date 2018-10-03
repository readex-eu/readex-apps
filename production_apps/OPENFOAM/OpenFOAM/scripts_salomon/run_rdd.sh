#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qexp
#PBS -N FOAM_RDD
#PBS -l select=1:ncpus=24:mpiprocs=24:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=01:00:00
#PBS -m be

APP='mpirun -n 24 simpleFoam -parallel'
PHASE_REG_NAME="iteration"

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
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
#readex-dyn-detect -p $PHASE_REG_NAME -t 0.1  scorep-*/profile.cubex
echo
echo "running readex-dyn-detect done"
