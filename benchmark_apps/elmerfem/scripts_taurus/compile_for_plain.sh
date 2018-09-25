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

SCRIPT_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
. ${SCRIPT_DIR}/readex_env/set_env_plain.source
. ${SCRIPT_DIR}/../paths.source

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
make install && echo "Elmer was successfully built."

cd ${ELMER_ROOT}
rm -f COMPILED_FOR_*
touch COMPILED_FOR_PLAIN
