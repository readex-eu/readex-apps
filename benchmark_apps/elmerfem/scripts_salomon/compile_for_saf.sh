#!/bin/sh

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_compile_saf
#PBS -l select=1:ncpus=24
#PBS -l walltime=02:00:00


if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
    LOC=$PBS_O_WORKDIR
else
    LOC=$(pwd)
fi

cd $LOC


# 1) Modules and variables
module load mkl
module load CMake/3.9.1

. readex_env/set_env_saf.source
. ../paths.source

cd ${ELMER_ROOT}




# 2) Compile
#export CC="scorep $READEX_CC"
#export CXX="scorep $READEX_CXX"
export FC="scorep $READEX_FC"
rm -rf install
rm -rf build
mkdir build
cd build
echo "Generate Make structure..."
cmake -DCMAKE_INCLUDE_CURRENT_DIR=ON -DCMAKE_INSTALL_PREFIX=../install .. && echo "Generating successful."
echo "Build Elmer."
make -j8 install && echo "Elmer was successfully built."
make -j4 install && echo "Elmer was successfully built."
make install && echo "Elmer was successfully built."

cd ${ELMER_ROOT}

rm -f COMPILED_FOR_*
touch COMPILED_FOR_SAF
