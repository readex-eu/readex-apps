#!/bin/sh

#SBATCH -t 05:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -A p_readex
#SBATCH --mem=24000
#SBATCH --reservation=READEX

# 1) Modules and variables
module load mkl/2017

SCRIPT_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
. ${SCRIPT_DIR}/readex_env/set_env_rdd.source
. ${SCRIPT_DIR}/../paths.source

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
