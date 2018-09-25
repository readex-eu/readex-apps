LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_meric.source
source $LOC/set_env_blasbench.source

make clean && make meric
cp blasbench blasbench_meric

