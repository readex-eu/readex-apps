LOC=$(pwd)
cd ..

module purge
source readex_env/set_env_ptf_rapl.source

export SCOREP_PREP="scorep --online-access --user --nocompiler --mpp=mpi --thread=omp --noopenmp"
export DSWITCH="-DUSE_SCOREP_MANUAL"

make clean
rm -rf lulesh2.0 lulesh2.0_ptf
make
mv lulesh2.0 lulesh2.0_ptf

