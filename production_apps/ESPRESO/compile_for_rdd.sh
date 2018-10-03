#!/bin/bash

#SBATCH -t 0-01:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -A p_readex
#SBATCH --mem-per-cpu=2500M
#SBATCH --reservation=READEX

source env/readex_env/set_env_rdd.source
source scripts_$READEX_MACHINE/environment.sh


# Check environment
module list
which $READEX_CXX
$READEX_CXX -v

# Compile ESPRESO
if [ $READEX_GNU ]
then
	FILTER_GCC="--instrument-filter=$(pwd)/scorep.filt"
	unset FILTER_ICC
	ls $(pwd)/scorep.filt
	if [ $? != 0 ]
	then
		echo "Compilation abort - missing scorep.filt"
		exit 1
	fi
else
	FILTER_ICC="-tcollect-filter=$(pwd)/scorep_icc.filt"
	unset FILTER_GCC
	ls $(pwd)/scorep_icc.filt
	if [ $? != 0 ]
	then
		echo "Compilation abort - missing scorep_icc.filt"
		exit 1
	fi
fi

rm -rf .waf-*
cp build.config.readex build.config
sed -i 's:$MERIC_ROOT:'$MERIC_ROOT':' build.config
sed -i 's/$READEX_CXX/'$READEX_CXX'/' build.config
sed -i 's/$READEX_CC/'$READEX_CC'/' build.config
sed -i 's/$READEX_OMP_FLAG/'$READEX_OMP_FLAG'/' build.config
sed -i 's:$READEX_INSTRUMENTATION:-DUSE_SCOREP '$FILTER_ICC':' build.config
sed -i 's:$READEX_SCOREP_FLAGS:--mpp=mpi --thread=omp --noopenmp --nomemory --online-access --user '$FILTER_GCC':' build.config
sed -i 's:$READEX_ATP_INCLUDE: :' build.config
sed -i 's:$READEX_ATP_LIB: :' build.config
sed -i 's/$READEX_ATP/ /' build.config

./waf distclean
./waf configure 
#./waf install 
./waf install -j 1 

mv espreso espreso_rdd
