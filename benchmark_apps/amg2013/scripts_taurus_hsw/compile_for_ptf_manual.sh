cd ..

module purge
source ./readex_env/set_env_ptf_hdeem.source
#source ./readex_env/set_env_ptf_rapl.source

export SCOREP_PREP="scorep --online-access --user --nocompiler --mpp=mpi --thread=none"

export DSWITCH="-DUSE_SCOREP_MANUAL"

make veryclean
make -j 12

cp ./test/amg2013 ./test/amg2013_ptf

