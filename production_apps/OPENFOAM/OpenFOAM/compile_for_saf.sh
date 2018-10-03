#SBATCH -t 00:15:00
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

source readex_env/set_env_saf.source
source scripts_$READEX_MACHINE/environment.sh

cd $FM_DIR/OpenFOAM-v1612+/

#edit make rules - wmake/rules/linux64Gcc/c++
export FOAM_CC="scorep --mpp=mpi --thread=none --nomemory g++ -std=c++11 -m64 -I$MERIC_ROOT/include"

source $FM_DIR/OpenFOAM-v1612+/etc/bashrc
export WM_NCOMPPROCS=24
export BOOST_ARCH_PATH=$BOOST_ROOT

cd applications/solvers/incompressible/simpleFoam/
wclean
wmake
