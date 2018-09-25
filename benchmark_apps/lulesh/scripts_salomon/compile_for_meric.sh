LOC=$(pwd)
cd ..

module purge
source readex_env/set_env_meric.source 

export DSWITCH="-DUSE_MERIC"

make clean
rm -rf lulesh2.0 lulesh2.0_meric
make
mv lulesh2.0 lulesh2.0_meric
