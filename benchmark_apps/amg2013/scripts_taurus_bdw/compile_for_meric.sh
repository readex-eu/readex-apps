cd ..

module purge
source ./readex_env/set_env_meric.source

export SCOREP_PREP=""
export DSWITCH="-DUSE_MERIC"

make veryclean
make -j 12

cp ./test/amg2013 ./test/amg2013_meric

