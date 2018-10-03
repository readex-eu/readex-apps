#!/bin/sh

#SBATCH -t 00:20:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -A p_readex
#SBATCH --mem-per-cpu=2500M

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
	export FM_DIR=$PBS_O_WORKDIR
else
	export FM_DIR=$(pwd)
fi
cd $FM_DIR

tar -zxvf OpenFOAM-v1612+.tgz
tar -zxvf ThirdParty-v1612+.tgz

source readex_env/set_env_saf.source
source scripts_$READEX_MACHINE/environment.sh

cd $FM_DIR/OpenFOAM-v1612+/
patch -p1 < ../patch.txt 
################################################################################
#PATCH:
# * applications/solvers/incompressible/simpleFoam/
#   	READEX instrumentration
# * wmake/rules/linux64Gcc/c++
#   	compilation flags
#   		export FOAM_CC=$(FOAM_CC)
#   		c++FLAGS    = $(GFLAGS) $(c++WARN) $(c++OPT) $(c++DBUG) $(ptFLAGS) $(LIB_HEADER_DIRS) -fPIC -fuse-ld=bfd
# * etc/bashrc
#   	FOAM_INST_DIR
################################################################################
# following sed commands should be edited if different settings is required
#sed -i 's/WM_COMPILER=Gcc/WM_COMPILER=Gcc/' OpenFOAM-v1612+/etc/bashrc
#sed -i 's/WM_MPLIB=SYSTEMOPENMPI/WM_MPLIB=SYSTEMOPENMPI/' OpenFOAM-v1612+/etc/bashrc

export FOAM_CC="g++ -std=c++11 -m64 -I$MERIC_ROOT/include"
source $FM_DIR/OpenFOAM-v1612+/etc/bashrc
export WM_NCOMPPROCS=24

export BOOST_ARCH_PATH=$BOOST_ROOT

# plain compilation of the OpenFOAM
./Allwmake

# problem decomposition
if [ $READEX_MACHINE == "taurus_bdw" ]
then
	cd $FM_DIR/../motorBike28/
else
	cd $FM_DIR/../motorBike24/
fi
chmod +x ./All*
./Allclean
./Allrun

