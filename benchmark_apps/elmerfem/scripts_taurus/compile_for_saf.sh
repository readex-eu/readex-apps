#!/bin/sh

#SBATCH -t 02:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -A p_readex
#SBATCH --mem=24000
#SBATCH --reservation=READEX

#Installation on Taurus:
# 1) Modules and variables
module load mkl/2017

SCRIPT_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
. ${SCRIPT_DIR}/readex_env/set_env_saf.source
. ${SCRIPT_DIR}/../paths.source

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
echo "Building Elmer..."
make -j8 install && echo "Elmer was successfully built."
make -j4 install && echo "Elmer was successfully built."
make install && echo "Elmer was successfully built."

cd ${ELMER_ROOT}

rm -f COMPILED_FOR_*
touch COMPILED_FOR_SAF
