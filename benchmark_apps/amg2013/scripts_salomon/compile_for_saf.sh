cd ..

module purge
source ./readex_env/set_env_saf.source

export SCOREP_PREP="scorep"
export DSWITCH=""

make veryclean
make -j 12 

cp ./test/amg2013 ./test/amg2013_saf

