module purge
source ./readex_env/set_env_rdd.source
source ./set_env_bem4i.source

export SCOREP_PREP="scorep --online-access --user --thread=omp --noopenmp --instrument-filter=scorep.filt"
export DSWITCH="-DUSE_SCOREP"

make CONF='release_taurus' clean
make -j 12 CONF='release_taurus' build

cp ./dist/release_taurus/bem4i ./dist/release_taurus/bem4i_rdd 
