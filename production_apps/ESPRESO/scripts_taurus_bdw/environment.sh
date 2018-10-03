#!/bin/sh

if [ $READEX_INTEL ]
then
	module load METIS/5.1.0-GCCcore-6.4.0
	module load CMake/3.11.4-GCCcore-6.4.0
	module load zlib/1.2.11
else
	module load METIS/5.1.0-GCCcore-7.3.0
	module load CMake/3.11.4-GCCcore-7.3.0
fi

source /projects/p_readex/it4i/mkl_2018_3/mkl/bin/mklvars.sh intel64

if [ "$READEX_GNU" == "1" ]; then
  export OMP_PLACES=cores
  export OMP_PROC_BIND=close
else
  export KMP_AFFINITY='compact,granularity=core'
fi

# mkl modules not add path to env - adding is done manualy
export CPATH=$MKL_INC:$CPATH
export LIBRARY_PATH=$MKL_LIB:$METIS_LIB:$LIBRARY_PATH
