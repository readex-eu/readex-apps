module purge
source ./readex_env/set_env_plain.source
source ./set_env_bem4i.source

export SCOREP_PREP=""
export DSWITCH=""

make CONF='release_taurus' clean
make -j 12 CONF='release_taurus' build

cp ./dist/release_taurus/bem4i ./dist/release_taurus/bem4i_plain 
