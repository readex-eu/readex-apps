LOC=$(pwd)
cd ..

module purge
source readex_env/set_env_saf.source

export SCOREP_PREP="scorep"
export DSWITCH="-DUSE_SCOREP"

make clean
rm -rf lulesh2.0 lulesh2.0_saf
make
mv lulesh2.0 lulesh2.0_saf

