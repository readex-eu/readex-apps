module purge
source ./readex_env/set_env_meric.source
source ./set_env_bem4i.source

export SCOREP_PREP=""
export DSWITCH="-DUSE_MERIC"

make CONF='release_taurus' clean
make -j 12 CONF='release_taurus' build

cp ./dist/release_taurus/bem4i ./dist/release_taurus/bem4i_meric 
