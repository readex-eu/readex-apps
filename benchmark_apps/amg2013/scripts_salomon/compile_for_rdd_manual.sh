cd ..

module purge
source ./readex_env/set_env_rdd.source

export SCOREP_PREP="scorep --online-access --user --nocompiler --mpp=mpi --thread=omp --noopenmp"
export DSWITCH="-DUSE_SCOREP_MANUAL"

make veryclean
make -j 12 

cp ./test/amg2013 ./test/amg2013_rdd

