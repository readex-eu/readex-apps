#!/bin/sh

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH -A p_readex
#SBATCH --mem=60000
#SBATCH -p haswell
#SBATCH --reservation=READEX


# 1) Modules and variables
module load mkl/2017
module load papi/5.5.1 

SCRIPT_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
. ${SCRIPT_DIR}/readex_env/set_env_meric.source
. ${SCRIPT_DIR}/../paths.source

cd ${ELMER_ROOT}


# 2) Compile
#export CXX="icc" 
#export FC="ifort" 
export FC="$READEX_FC"
export CXX="$READEX_CXX"

${FC} -c ${MERIC_ROOT}/include/meric_mod.F90

rm -rf install
rm -rf build
mkdir build
cd build
echo "Generate Make structure..."
cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE -DCOMPILE_FOR_MERIC:BOOL=TRUE -DCMAKE_INCLUDE_CURRENT_DIR=ON -DCMAKE_INSTALL_PREFIX=../install .. && echo "Generating successful."
#cmake -DCMAKE_INCLUDE_CURRENT_DIR=ON -DCMAKE_INSTALL_PREFIX=../install .. && echo "Generating successful."
echo "Build Elmer."
make -j8 install && echo "Elmer was successfully built."
make -j4 install && echo "Elmer was successfully built."
make install && echo "Elmer was successfully built."

cd ${ELMER_ROOT}
rm -f COMPILED_FOR_*
touch COMPILED_FOR_MERIC

