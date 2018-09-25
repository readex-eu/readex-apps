LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_plain.source
source $LOC/set_env_blasbench.source

make clean && make
cp blasbench blasbench_plain

