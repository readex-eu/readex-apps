#!/bin/sh

#SBATCH -t 30:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -A p_readex
#SBATCH --mem=62000

#Installation on Taurus:
# 1) Modules and variables
. readex_env/set_env_meric.source
. scripts_$READEX_MACHINE/environment.sh

# 2) Compile
cp CMakeLists-MERIC.txt CMakeLists.txt
export CXX=$READEX_CXX
rm -rf build
mkdir build
cd build
cmake ..
make


# Note: Chose MERIC or SCOREP by adding -DUSE_MERIC or -DUSE_SCOREP, respectively, 
#       in CMakeLists.txt, variable CMAKE_CXX_FLAGS

# Note: READEX kernels are in src/kripke.cpp (Main region) and
#       src/Kripke/Sweep_Solver.cpp (Compute regions)
 
