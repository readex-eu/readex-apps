#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_compile_plain
#PBS -l select=1:ncpus=24

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
    LOC=$PBS_O_WORKDIR
else
    LOC=$(pwd)
fi

cd $LOC


# 1) Modules and variables
module load mkl
module load CMake/3.9.1

. readex_env/set_env_plain.source
. ../paths.source

cd ${ELMER_ROOT}


# 2) Compile
export CXX=$READEX_CXX
export FC=$READEX_FC
rm -rf build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install ..
make -j8 install
make -j4 install
make install && echo "Elmer successfully built."

rm -f COMPILED_FOR_*
touch COMPILED_FOR_PLAIN
