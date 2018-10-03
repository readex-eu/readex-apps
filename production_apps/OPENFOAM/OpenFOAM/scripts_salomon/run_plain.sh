#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qexp
#PBS -N FOAM
#PBS -l select=1:ncpus=24:mpiprocs=24:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=01:00:00
#PBS -m be

################################################################################
if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
export FM_DIR=$FM_DIR/..

cd $FM_DIR
source readex_env/set_env_saf.source
source scripts_$READEX_MACHINE/environment.sh
cd $FM_DIR/OpenFOAM-v1612+/
source $FM_DIR/OpenFOAM-v1612+/etc/bashrc
################################################################################

#cd tutorials/incompressible/simpleFoam/motorBike24/
#mpirun -n 24 simpleFoam -parallel
cd ../../motorBike24/
./Allclean
./Allrun
