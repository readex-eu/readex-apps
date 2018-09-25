#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_compile_rdd
#PBS -l select=1:ncpus=24
#PBS -l walltime=01:00:00


if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
    LOC=$PBS_O_WORKDIR
else
    LOC=$(pwd)
fi

cd $LOC


# 1) Modules and variables
module load mkl
module load CMake/3.9.1

. readex_env/set_env_rdd.source
. ../paths.source


cd ${ELMER_ROOT}


# 2) Compile
if [ $READEX_INTEL ]
then
        FILTER_ICC="-tcollect-filter=$(pwd)/scorep_icc.filt"
    unset FILTER_GCC
else
        FILTER_GCC="--instrument-filter=$(pwd)/scorep.filt"
    unset FILTER_ICC
fi

export FC="scorep --online-access --user --mpp=mpi --thread=none --nomemory $FILTER_GCC $READEX_FC $FILTER_ICC"
rm -rf install
rm -rf build
mkdir build
cd build
cmake -DCOMPILE_FOR_RDD:BOOL=TRUE -DCMAKE_INSTALL_PREFIX=../install ..
make -j8 install
make -j4 install
make install && echo "Elmer build complete."

cd ${ELMER_ROOT}
rm -f COMPILED_FOR_*
touch COMPILED_FOR_RDD
