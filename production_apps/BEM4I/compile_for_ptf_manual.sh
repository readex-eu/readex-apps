module purge
source ./readex_env/set_env_ptf_hdeem.source
#source ./readex_env/set_env_ptf_rapl.source
source ./set_env_bem4i.source

export SCOREP_PREP="scorep --online-access --user --nocompiler --thread=omp --noopenmp"
export DSWITCH="-DUSE_SCOREP_MANUAL"

make CONF='release_taurus' clean
make -j 12 CONF='release_taurus' build

cp ./dist/release_taurus/bem4i ./dist/release_taurus/bem4i_ptf 
