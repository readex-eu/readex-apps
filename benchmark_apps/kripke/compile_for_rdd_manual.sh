#!/bin/sh

#SBATCH -t 30:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -A p_readex
#SBATCH --mem=62000

#Installation on Taurus:
# 1) Modules and variables
. readex_env/set_env_rdd.source
. scripts_$READEX_MACHINE/environment.sh

# 2) Compile
cp CMakeLists-SCOREP-MANUAL.txt CMakeLists.txt
export CXX="scorep --online-access --user --mpp=mpi --thread=none --nomemory --nocompiler $READEX_CXX"
rm -rf build
mkdir build
cd build
cmake ..
make

cp ../scripts_$READEX_MACHINE/run_rdd.sh .

# 3) Test Run 

#srun -n 24 ./kripke --procs 2,2,6 --niter 20 --nest GZD --zones 4,4,4 --groups 200
#srun -n 24 ./kripke --procs 2,2,6 --niter 10 --nest GZD --zones 32,32,32 --legendre 8 --dset 32


# Note: Chose MERIC or SCOREP by adding -DUSE_MERIC or -DUSE_SCOREP, respectively, 
#       in CMakeLists.txt, variable CMAKE_CXX_FLAGS

# Note: READEX kernels are in src/kripke.cpp (Main region) and
#       src/Kripke/Sweep_Solver.cpp (Compute regions)
 
