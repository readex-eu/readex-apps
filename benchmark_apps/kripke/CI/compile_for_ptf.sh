#!/bin/sh

#SBATCH -t 30:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -A p_readex
#SBATCH --mem=62000

#Installation on Taurus:
# 1) Modules and variables
cd ..
. readex_env/set_env_ptf_hdeem.source
. $(pwd)/environment.sh

# 2) Compile
cp CMakeLists-SCOREP-PHASE.txt CMakeLists.txt
if [ $READEX_INTEL ]
then
        FILTER_ICC="-tcollect-filter=$(pwd)/RESULTS/scorep_icc.filt"
else
        FILTER_GCC="--instrument-filter=$(pwd)/RESULTS/scorep.filt"
fi
export CXX="scorep --online-access --user --mpp=mpi --thread=none --nomemory $FILTER_GCC $READEX_CXX $FILTER_ICC"
rm -rf build
mkdir build
cd build
cmake ..
make

cp ../scripts/run_ptf.sh .
cp ../scripts/run_rrl.sh .


# 3) Test Run 

#srun -n 24 ./kripke --procs 2,2,6 --niter 20 --nest GZD --zones 4,4,4 --groups 200
#srun -n 24 ./kripke --procs 2,2,6 --niter 10 --nest GZD --zones 32,32,32 --legendre 8 --dset 32


# Note: Chose MERIC or SCOREP by adding -DUSE_MERIC or -DUSE_SCOREP, respectively, 
#       in CMakeLists.txt, variable CMAKE_CXX_FLAGS

# Note: READEX kernels are in src/kripke.cpp (Main region) and
#       src/Kripke/Sweep_Solver.cpp (Compute regions)
 
