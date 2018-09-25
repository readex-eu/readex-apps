module purge
source ./readex_env/set_env_saf.source
source ./set_env_bem4i.source

export SCOREP_PREP="scorep"
export DSWITCH=""

make CONF='release_taurus' clean
make -j 12 CONF='release_taurus' build

cp ./dist/release_taurus/bem4i ./dist/release_taurus/bem4i_saf 
